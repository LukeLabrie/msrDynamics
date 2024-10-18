from jitcdde import jitcdde, y, t, jitcdde_input, input
import numpy as np
import chspy
import sympy as sp
import matplotlib.pyplot as plt

class System:
    """
    System class for thermal hydraulic system.
    
    This class represents the entire thermal hydraulic system, consisting of multiple nodes.

    Attributes:
        nodes (dict): Dictionary of nodes in the system.
        dydt (list): List of symbolic expressions representing the rate of change for each node.
        y0 (list): List of initial conditions for the system.
        inputs (None or CubicHermiteSpline): Input spline for the system.
        input_funcs (list): List of input functions.
        integrator (jitcdde or jitcdde_input): JiTCDDE integrator.
        custom_past (None or CubicHermiteSpline): Custom past spline for the integrator.
        trip_conditions (None or list): List of trip conditions.
        trip_info (None or dict): Dictionary containing trip information.
    """

    def __init__(self, nodes=None, dydt=None, y0=None) -> None:
        """
        Initialize the System.
        
        Parameters:
            nodes (dict, optional): Dictionary of nodes in the system. Default is None.
            dydt (list, optional): List of symbolic expressions for the rate of change. Default is None.
            y0 (list, optional): List of initial conditions. Default is None.
        """
        self.n_nodes = 0
        self.n_inputs = 0

        # avoids using mutable objects as default arguments, which causes a memory leak
        if not nodes:
            self.nodes = {}
        else:
            self.nodes = nodes
        if not dydt:
            self._dydt = []
        else:
            self._dydt = dydt
        if not y0:
            self.y0 = []
        else:
            self.y0 = y0

        self.inputs = None
        self.input_funcs = None 
        self.integrator = None
        self.custom_past = None
        self.trip_conditions = None
        self.trip_info = {
            'time': float('inf'),
            'tripped': False
            }
        self.input = None
        self.pid_controllers = None
        self.input_func_names = None

    @property
    def dydt(self):
        """list: List of symbolic expressions representing the rate of change for each node."""
        return [node.dydt for node in self.nodes.values()]

    @dydt.setter
    def dydt(self, dydt):
        self._dydt = dydt

    def _get_full_input(self, times):
        return [f(times) for f in self.input_funcs]

    def _check_trip(self, time, state, diff):
        """
        Check if any trip condition is met.

        Parameters:
            time (float): Current time.
            state (list): Current state values.
            diff (list): Current derivative values.

        Returns:
            tuple: Index of tripped condition and trip limit if tripped, None otherwise.
        """
        for idx, trip_obj in enumerate(self.trip_conditions):
            s = None
            if trip_obj.check_after:
                if time < trip_obj.check_after:
                    return None
            bounds_s = trip_obj.bounds
            min_s = bounds_s[0]
            max_s = bounds_s[1]

            # if delayed, interpolate past to check trip condition
            if trip_obj.delay:
                if (len(self.trip_info['state']) > 1) and (time >= trip_obj.delay):
                    chs = self.trip_info['state']
                    idx_interp = chs.last_index_before(time-trip_obj.delay)
                    a1 = chs[idx_interp]
                    a2 = chs[idx_interp+1]
                    if trip_obj.trip_type == 'state':
                        s = chspy.interpolate(time - trip_obj.delay, idx, (a1,a2))
                    elif trip_obj.trip_type == 'diff':
                        s = chspy.interpolate_diff(time - trip_obj.delay, idx, (a1,a2))
                    else:
                        raise ValueError('''Invalid trip type. Currently supported 
                                            are 'state' and 'diff'.''')
            else:
                if trip_obj.trip_type == 'state':
                    s = state[idx]
                elif trip_obj.trip_type == 'diff':
                    s = diff[idx]
                else:
                    raise ValueError('''Invalid trip type. Currently supported 
                                       are 'state' and 'diff'.''')
            if s:
                if (s < min_s):
                    return (idx, min_s)
                elif (s > max_s):
                    return (idx, max_s)
        return None

    def add_input(self, input_func, times, max_anchors = 100, input_tol = 5, name: str = None):
        """
        Add an input function of time to the integrator.

        Parameters:
            input_func (callable): Input function.

        Returns:
            input: JiTCDDE input object.
        """
        # store input functions
        if self.input_funcs:
            self.input_funcs.append(input_func)
            self.input_func_names.append(name)
        else:
            self.input_funcs = [input_func]
            self.input_func_names = [name]

        # set up input spline
        spline = chspy.CubicHermiteSpline(n=1)
        spline.from_function(input_func, 
                                times_of_interest=times, 
                                max_anchors=max_anchors, 
                                tol=input_tol)
        if self.input:
            self.input = chspy.join(self.input,spline)
        else:
            self.input = spline

        # return jitcdde input object so it can be used for sympy expressions
        input_obj = input(self.n_inputs)
        self.n_inputs += 1

        return input_obj

    def set_custom_past(self, past: chspy.CubicHermiteSpline, t_truncate=None):
        """
        Set custom past data for the integrator.

        Parameters:
            past (chspy.CubicHermiteSpline): Custom past data.
            t_truncate (float, optional): Time to truncate the past data. Default is None.
        """
        if t_truncate:
            past.truncate(t_truncate)
        self.custom_past = past

    def finalize(self, T, sdd=False, md=1e10):
        """
        Instantiate and store JiTCDDE integrator.

        Parameters:
            T (np.ndarray): Time array.
            sdd (bool): State-dependent delays.
            md (float): Maximum delay.
            max_anchors (int): Maximum number of anchors.
            input_tol (int): Tolerance for input spline.
        """
        if self.integrator: 
            self.integrator = None
        # set up system matrix
        self.dydt = [n.dydt for n in self.nodes.values()]
        if self.input:
            # instantiate integrator
            if self.pid_controllers:
                DDE = jitcdde_input(self.dydt, 
                                    self.input, 
                                    max_delay=md, 
                                    callback_functions=[ (pc.sym_func, pc.func, pc.n_args) for pc in self.pid_controllers])
            else:
                DDE = jitcdde_input(self.dydt, self.input, max_delay=md)

            # set initial conditions
            if not self.custom_past:
                self.y0 = [n.y0 for n in self.nodes.values()]
                DDE.constant_past(self.y0)
            else:
                DDE.purge_past()
                # shift past time to end at t = 0
                t_last = self.custom_past[-1].time
                for a in self.custom_past:
                    a.time -= t_last
                DDE.add_past_points(self.custom_past)
            DDE.adjust_diff()

        else:
            if self.pid_controllers:
                DDE = jitcdde(self.dydt, max_delay=md, callback_functions=[ (pc.sym_func, pc.func, pc.n_args) for pc in self.pid_controllers])
            else:
                DDE = jitcdde(self.dydt, max_delay=md)
            # set initial conditions
            if not self.custom_past:
                self.y0 = [n.y0 for n in self.nodes.values()]
                DDE.constant_past(self.y0)
            else:
                # shift past time to end at t = 0
                t_last = self.custom_past[-1].time
                for a in self.custom_past:
                    a.time -= t_last
                DDE.add_past_points(self.custom_past)

        # max delay needs to be provided in the case of state-dependent delays
        if sdd:
            DDE.max_delay = md
        self.integrator = DDE

    def add_nodes(self, new_nodes: list):
        """
        Add nodes to the system.

        Parameters:
            new_nodes (list): List of nodes to be added.
        """
        index = self.n_nodes
        for n in new_nodes: 
            n.y = lambda tau=None, i=index: y(i, tau) if tau is not None else y(i)
            n.index = index
            # add node to nodes dictionary by name, if available, and by index if not
            if n.name:
                assert n.name not in self.nodes, f"Nodes cannot share the same name. {n.name} is a duplicate."
                self.nodes[n.name] = n
            else:
                self.nodes[index] = n
            self.n_nodes += 1
            index += 1

    def get_state_by_index(self, i: int, j: int):
        """
        Get the state of a node by its index.

        Parameters:
            i (int): Index of the node.
            j (int): Time index.

        Returns:
            tuple: State value and derivative at the specified index.
        """
        val = self.integrator.get_state()[j][1][i]
        deriv = self.integrator.get_state()[j][2][i]
        return (val, deriv)

    def solve(self, 
              T, 
              sdd=False, 
              max_delay=1e10, 
              populate_nodes=False, 
              abs_tol=1e-10, 
              rel_tol=1e-05, 
              min_step = 1e-10, 
              max_step = 10.0):
        """
        Solve the system and return the solution matrix.

        Parameters:
            t0 (float): Initial time.
            tf (float): Final time.
            incr (float): Time increment.
            sdd (bool): State-dependent delays.
            max_delay (float): Maximum delay.
            max_anchors (int): Maximum number of anchors.
            input_tol (int): Tolerance for input spline.
            populate_nodes (bool): Populate node objects with solutions.
            abs_tol (float): Absolute tolerance for integration.
            rel_tol (float): Relative tolerance for integration.

        Returns:
            np.ndarray: Solution matrix.
        """
        if len(self.nodes) == 0:
            raise ValueError('No nodes have been added to the system')
        
        # clear data
        for n in self.nodes.values():
            if (n.y_out.any()):
                n.y_out = []

        # set integrator 
        print("finalizing integrator...")
        self.finalize(T, sdd, max_delay)
        self.integrator.set_integration_parameters(atol=abs_tol, rtol=rel_tol, min_step = min_step, max_step = max_step)

        # solution 
        y = []

        print("integrating...")
        if self.trip_conditions:

            # integrate with trip conditions
            for t_x in T:
                # extract state and derivs for trip check 
                y.append(np.array(self.integrator.integrate(t_x)))
                idxs = [c.idx for c in self.trip_conditions]
                states = y[-1][idxs]

                # derivative is only estimated after the first step 
                if len(y) > 1:
                    derivs = ((y[-1]-y[-2])/(T[1]-T[0]))[idxs]
                else:
                    derivs = [0.0]*len(states)
                
                if 'state' in self.trip_info:
                    self.trip_info['state'].extend([chspy.Anchor(t_x, states, derivs)])
                else:
                    self.trip_info['state'] = chspy.CubicHermiteSpline(n=len(self.trip_conditions), 
                                                          anchors=[chspy.Anchor(t_x, states, derivs)])

                # check if system has tripped
                tripped = self._check_trip(t_x, states, derivs)
                if tripped:
                    # get trip condition object
                    trip_obj = self.trip_conditions[tripped[0]]
                    print(f'idx {tripped[0]} tripped after integration to t = {t_x:3f} with a value of {tripped[1]}')

                    # store trip info 
                    self.trip_info['tripped'] = True
                    self.trip_info['idx'] = trip_obj.idx
                    self.trip_info['limit'] = tripped[1]
                    self.trip_info['type'] = trip_obj.trip_type

                    # get system spline
                    print('getting state...')
                    state = self.integrator.get_state()

                    # calculate exact trip time using splines 
                    print('computing trip time within interval...')
                    trip_sol = []
                    start = trip_obj.check_after if trip_obj.check_after is not None else state[0].time
                    solve_diff = True if self.trip_info['type'] == 'diff' else False
                    trip_sol = state.solve(self.trip_info['idx'],
                                            self.trip_info['limit'],
                                            beginning=start,
                                            solve_derivative = solve_diff)
                    if trip_obj.delay:
                        self.trip_info['time'] = trip_sol[0][0] + trip_obj.delay
                    else:
                        self.trip_info['time'] = trip_sol[0][0] 
                    print(f"tripped at t = {self.trip_info['time']:.3f}")
                    print(f"state idx: {tripped[0]}")
                    print(f"limit: {tripped[1]}")
                    break
        else:
            for t_x in T:
                y.append(self.integrator.integrate(t_x))

        # populate node objects with solutions 
        if populate_nodes:
            print('populating nodes objects solution vectors...')
            for s in enumerate(self.nodes.values()):
                s[1].y_out = np.array([state[s[0]] for state in y])

        return np.array(y)

    def plot_input(self, index):
        """
        Plot the input function for the given index.

        Parameters:
            index (int): Index of the input function.

        Returns:
            matplotlib.axes.Axes: Axes object with the plot of the selected input.
        """
        if self.input is None:
            raise ValueError('No input has been defined')
        else:
            _, ax = plt.subplots()
            times = []
            vals = []
            for i in self.input:
                times.append(i.time)
                vals.append(i.state[index])
            ax.plot(times, vals)
            return ax


class Node:
    """
    Node class for thermal hydraulic system.
    
    This class represents a node in an MSR system, encapsulating physical properties and dynamics.

    Attributes:
        name (str): The name of the node.
        m (float): Mass of the node in kilograms (kg).
        scp (float): Specific heat capacity of the node in joules per kilogram kelvin (J/(kg*K)).
        W (float): Mass flow rate through the node in kilograms per second (kg/s).
        y0 (float): Initial temperature of the node in kelvin (K).
        dTdt_advective (float): Symbolic expression for the rate of temperature change due to advective heat flow in kelvin per second (K/s).
        dTdt_internal (float): Symbolic expression for the rate of temperature change due to internal heat generation in kelvin per second (K/s).
        dTdt_convective (float): Symbolic expression for the rate of temperature change due to convective heat flow in kelvin per second (K/s).
        dndt (float): Symbolic expression for the rate of change of neutron population.
        dcdt (float): Symbolic expression for the rate of change of precursor concentration.
        drdt (float): Symbolic expression for the rate of change of reactivity.
        y (None or callable): JiTCDDE state variable object, to be assigned by the System class.
        y_out (np.ndarray): Solution data, to be populated by the System class.
    """

    def __init__(self, name: str = None, m: float = 0.0, scp: float = 0.0, W: float = 0.0, y0: float = 0.0) -> None:
        """
        Initialize the node.
        
        Parameters:
            name (str): The name of the node.
            m (float): Mass of the node in kilograms (kg).
            scp (float): Specific heat capacity of the node in joules per kilogram kelvin (J/(kg*K)).
            W (float): Mass flow rate through the node in kilograms per second (kg/s).
            y0 (float): Initial temperature of the node in kelvin (K).
        """
        self.name = name            # node name
        self.m = m                  # mass (kg)
        self.scp = scp              # specific heat capacity (J/(kg*°K))
        self.W = W                  # mass flow rate (kg/s)
        self.y0 = y0                # initial temperature (°K)
        self.dTdt_advective = 0.0   # sym. expression for advective heat flow (°K/s)
        self.dTdt_internal = 0.0    # sym. expression for internal heat generation (°K/s)
        self.dTdt_convective = 0.0  # sym. expression for convective heat flow (°K/s)
        self.dndt = 0.0             # sym. expression for dn/dt (n = neutron population)
        self.dcdt = 0.0             # sym. expression for dc/dt (c = precursor concentration)
        self.drdt = 0.0             # sym. expression for dr/dt (r = reactivity)
        self.dndt_decay = 0.0       # syms expression for dn_d/dt (n_d = fractional decay generation)
        self._dydt = 0.0            # sym. expression for user-defined dynamics 
        self.y = None               # JiTCDDE state variable object, to be assigned by System
        self.index = None           # JiTCDDE state variable index, to be assigned by System
        self.y_out = np.array([])   # solution data, to be populated by System
        self.y_rhs = np.array([])   # solution data, to be populated by System
        self.trip_conditions = None

    @property
    def dydt(self):
        """float: Symbolic expression for the rate of change of state variables."""
        return self.dTdt_advective + self.dTdt_internal + \
               self.dTdt_convective + self.dndt + self.dcdt + self.drdt + \
               self._dydt + self.dndt_decay

    @dydt.setter
    def dydt(self, custom_dydt):
        self._dydt = custom_dydt

    def set_dTdt_advective(self, source):
        """
        Set the rate of temperature change due to advective heat transfer.

        Parameters:
            source (Node or float): Source node or constant temperature.
        """
        # check not to mix thermal and point-kinetic equations
        if self.dndt or self.dcdt or self.drdt:
            raise ValueError('''This node has already been assigned 
                             point-kinetic dynamics''')
        
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dTdt_advective = 0.0
            self.dTdt_advective = (source - self.y()) * self.W / self.m
        else:
            raise ValueError('''Nodes need to be added to a System() object 
                             before setting dynamics.''')

    def set_dTdt_internal(self, source: list, k: list):
        """
        Set the rate of temperature change due to internal heat generation.

        Parameters:
            source (callable): Source of generation (state variable y(i)).
            k (float): Constant of proportionality.
        """
        # check not to mix thermal and point-kinetic equations
        if self.dndt or self.dcdt or self.drdt:
            raise ValueError('''This node has already been assigned 
                             point-kinetic dynamics''')
        
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dTdt_internal = 0.0
            for idx, s in enumerate(source):
                self.dTdt_internal += k[idx] * s / (self.m * self.scp)
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

    def set_dTdt_convective(self, source: list, hA: list):
        """
        Set the rate of temperature change due to convective heat transfer.

        Parameters:
            source (list): List of source state variables y(i).
            hA (list): List of convective heat transfer coefficients * wetted areas (MW/C).
        """
        # check not to mix thermal and point-kinetic equations
        if self.dndt or self.dcdt or self.drdt:
            raise ValueError('''This node has already been assigned 
                             point-kinetic dynamics''')
        if (self.m <= 0.0) or (self.scp <= 0.0):
            print(f'node mass: {self.m:.2f}')
            print(f'node specific heat capacity: {self.scp:.2f}')
            raise ValueError('''Invalid specific heat capacity or mass for node. To set convective heat transfer 
                                dynamics, both mass and specific heat capacity need to be greater than zero.''')
        
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dTdt_convective = 0.0
            
            # make sure there is an associated coefficient for each source
            if len(source) != len(hA):
                raise ValueError("Source array and coefficient array different lengths")
            for i in range(len(source)):
                self.dTdt_convective += hA[i] * (source[i] - self.y()) / (self.m * self.scp)
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

    def set_dndt(self, r: y, beta_eff: float, Lambda: float, lam: list, C: list):
        """
        Set the rate of change of neutron population using the point kinetics equation.

        Parameters:
            r (y): Reactivity, rho (state variable y(i)).
            beta_eff (float): Delayed neutron fraction.
            Lambda (float): Prompt neutron lifetime (s).
            lam (list): Decay constants for associated precursor groups.
            C (list): Precursor groups (list of state variables y(i)).
        """
        # check not to mix thermal and point-kinetic equations
        if self.dTdt_advective or \
           self.dTdt_internal or \
           self.dTdt_convective or \
           self.drdt or \
           self.dcdt:
            raise ValueError('''This node has already been assigned 
                             incompatible dynamics''')
        
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dndt = 0.0
            precursors = 0.0
            for g in enumerate(lam):
                precursors += g[1] * C[g[0]]
            self.dndt = (r - beta_eff) * self.y() / Lambda + precursors
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

    def set_dcdt(self, n: y, beta: float, Lambda: float, lam: float, flow: bool = False, t_c: float = 0.0, t_l: float = 0.0):
        """
        Set the rate of change of precursor concentration.

        Parameters:
            n (y): Fractional neutron concentration (state variable y(i)).
            beta (float): Fraction for precursor group.
            Lambda (float): Prompt neutron lifetime.
            lam (float): Decay constant for precursor group.
            flow (bool): Indicates if flow is considered in the model.
            t_c (float): Core transit time (seconds).
            t_l (float): Loop transit time (seconds).

        Raises:
            ValueError: If flow is True and either t_c or t_l are not set properly.
        """
        # check not to mix thermal and point-kinetic equations
        if self.dTdt_advective or \
           self.dTdt_internal or \
           self.dTdt_convective or \
           self.drdt or \
           self.dndt:
            raise ValueError('''This node has already been assigned 
                             incompatible dynamics''')
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dcdt = 0.0
            source = n * beta / Lambda
            decay = lam * self.y()
            if flow:
                outflow = self.y() / t_c
                inflow = self.y(t - t_l) * sp.exp(-lam * t_l) / t_c
                self.dcdt = source - decay - outflow + inflow
            else:
                self.dcdt = source - decay
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

    def set_drdt(self, sources: list, coeffs: list):
        """
        Set the rate of change of reactivity based on feedback.

        Parameters:
            sources (list): List of derivatives of feedback sources (dy(i)/dt).
            coeffs (list): List of respective feedback coefficients.
        """
        # check not to mix thermal and point-kinetic equations
        if self.dTdt_advective or \
           self.dTdt_internal or \
           self.dTdt_convective or \
           self.dndt or \
           self.dcdt:
            raise ValueError('''This node has already been assigned 
                             incompatible dynamics''')
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.drdt = 0.0
            fb = 0.0
            for s in enumerate(sources):
                fb += s[1] * coeffs[s[0]]
            self.drdt = fb
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")
        
    def set_dndt_decay(self, n: y, n0: float, rel_yield: float, lam: float):
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dndt_decay = 0.0
            self.dndt_decay += ((n/n0) * rel_yield - lam * self.y())  
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")
