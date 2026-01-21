from jitcdde import jitcdde, y, t, jitcdde_input, input, current_y
import numpy as np
import chspy
import sympy as sp
import matplotlib.pyplot as plt
from tqdm import tqdm
from symengine import Mul
import sys
import gc
from symengine import lambdify
import numpy as np
import symengine as se
import os

MAX_INT = sys.maxsize

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
        self.pid_loops = []
        self.input_func_names = None
        self.max_delay = None
        self.callback_functions = []

    @property
    def dydt(self):
        """list: List of symbolic expressions representing the rate of change for each node."""
        return [node.dydt for node in self.nodes.values()]

    @dydt.setter
    def dydt(self, dydt):
        self._dydt = dydt

    @property
    def y0(self):
        """list: List of intial conditions for each node."""
        return np.array([node.y0 for node in self.nodes.values()])
    
    @y0.setter
    def y0(self, y0):
        for idx, n in enumerate(self.nodes):
            self.nodes[n].y0 = y0[idx]

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
                    elif trip_obj.trip_type == 'diff_rel':
                        s = chspy.interpolate_diff(time - trip_obj.delay, idx, (a1,a2))/chspy.interpolate(time - trip_obj.delay, idx, (a1,a2))
                    else:
                        raise ValueError('''Invalid trip type. Currently supported 
                                            are 'state' and 'diff'.''')
            else:
                if trip_obj.trip_type == 'state':
                    s = state[idx]
                elif trip_obj.trip_type == 'diff':
                    s = diff[idx]
                elif trip_obj.trip_type == 'diff_rel':
                    s = diff[idx]/state[idx]
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

    def finalize(self, start_time):
        """
        Instantiate and store JiTCDDE integrator.

        Parameters:
            T (np.ndarray): Time array.
            sdd (bool): State-dependent delays.
            md (float): Maximum delay.
            max_anchors (int): Maximum number of anchors.
            input_tol (int): Tolerance for input spline.
        """
        # set up system matrix
        self.dydt = [n.dydt for n in self.nodes.values()]

        # clang results in faster integration
        os.environ["CC"] = "clang"

        # input uses different integrator object
        if self.input:
            DDE = jitcdde_input(self.dydt, self.input, max_delay = self.max_delay)
        else:
            DDE = jitcdde(self.dydt, max_delay = self.max_delay)

        # populate callback functions 
        if self.pid_loops:
            self.callback_functions.extend([ (pc.output_sym, pc.output_func, pc.n_args) for pc in self.pid_loops])
        DDE.callback_functions = self.callback_functions

        # set initial conditions
        if not self.custom_past:
            self.y0 = [n.y0 for n in self.nodes.values()]
            DDE.constant_past(self.y0, time = start_time)
        else:
            # shift past time to end at t = 0
            t_last = self.custom_past[-1].time
            for a in self.custom_past:
                a.time -= t_last
            DDE.add_past_points(self.custom_past)

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
            n.in_system = True

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
    
    def get_node_by_index(self, idx):

        for _, n in self.nodes.items():
            if n.index == idx:
                return n
        raise ValueError(f'Node with index {idx} not found')
    
    def _prepare_integrator(self, 
                            start_time=0.0,
                            abs_tol=1e-10, 
                            rel_tol=1e-05, 
                            min_step = 1e-10, 
                            max_step = 10.0,
                            ):
        
        if len(self.nodes) == 0:
            raise ValueError('No nodes have been added to the system')
        
        # clear data
        for n in self.nodes.values():
            if (n.y_out.any()):
                n.y_out = []

        # set integrator 
        print("finalizing integrator...")
        self.finalize(start_time)
        self.integrator.set_integration_parameters(atol=abs_tol, rtol=rel_tol, min_step = min_step, max_step = max_step)
    
    def solve(self, 
              T, 
              max_delay=1e10, 
              populate_nodes=False, 
              abs_tol=1e-10, 
              rel_tol=1e-05, 
              min_step = 1e-10, 
              max_step = 10.0,
              md_step = 1e-3,
              ):
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
        self.max_delay = max_delay
        self._prepare_integrator(T[0], abs_tol, rel_tol, min_step, max_step)

        print("integrating...")
        # solve
        if self.trip_conditions:
            y = self._solve_with_trip_conditions(T, md_step)
        else:
            y = self._solve_default(T, md_step)

        # populate node objects with solutions, off by default, as it can cause
        # memory blowup/leakage when running many models 
        if populate_nodes:
            self._populate_nodes(y)   

        return np.array(y)

    def _solve_default(self, times, md_step):
        y = []
        with tqdm(total=len(times), desc="Integration progress") as pbar:
            for t_x in times[times<=(self.max_delay+times[0])]:
                y.append(self.integrator.integrate_blindly(t_x, step = md_step))
                pbar.update(1)
            for t_x in times[times>(self.max_delay+times[0])]:
                y.append(self.integrator.integrate(t_x))
                pbar.update(1)
        return np.array(y)

    def _solve_with_trip_conditions(self, times, md_step):
        y = []
        # integrate with trip conditions
        with tqdm(total=len(times), desc="Integration progress") as pbar:
            for t_idx, t_x in enumerate(times):
                # extract state and derivs for trip check 
                if (t_x <= (self.max_delay+times[0])):
                    y.append(np.array(self.integrator.integrate_blindly(t_x, step = md_step)))
                else:
                    y.append(np.array(self.integrator.integrate(t_x)))
                idxs = [c.idx for c in self.trip_conditions]
                states = y[-1][idxs]

                # derivative is only estimated after the first step 
                if len(y) > 1:
                    derivs = (((y[-1]-y[-2])/(times[t_idx]-times[t_idx-1])))[idxs]
                else:
                    derivs = [0.0]*len(states)
                
                if 'state' in self.trip_info:
                    self.trip_info['state'].extend([chspy.Anchor(t_x, states, derivs)])
                    self.trip_info['state'].forget(self.max_delay)
                else:
                    self.trip_info['state'] = chspy.CubicHermiteSpline(n=len(self.trip_conditions), 
                                                        anchors=[chspy.Anchor(t_x, states, derivs)])

                # check if system has tripped
                tripped = self._check_trip(t_x, states, derivs)
                pbar.update(1)
                if tripped:
                    # get trip condition object
                    trip_obj = self.trip_conditions[tripped[0]]
                    print(f'{trip_obj.name} tripped after integration to t = {t_x:3f} with a value of {tripped[1]}')

                    # store trip info 
                    self.trip_info['tripped'] = True
                    self.trip_info['idx'] = trip_obj.idx
                    self.trip_info['limit'] = tripped[1]
                    self.trip_info['type'] = trip_obj.trip_type
                    self.trip_info['t_last'] = t_x
                    self.trip_info['name'] = trip_obj.name

                    # get system spline
                    print('getting state...')
                    state = self.integrator.get_state()

                    # calculate exact trip time using splines 
                    print('computing trip time within interval...')
                    trip_sol = []

                    # bounds for inteprolation

                    if self.trip_info['type'] == 'diff_rel':
                        # set up new spline for fractional derivative and interpolate
                        times_dr = np.array([s.time for s in state])
                        states_dr = np.array([s.diff[self.trip_info['idx']]/s.state[self.trip_info['idx']] for s in state])
                        deriv_dr = np.diff(states_dr)/np.diff(times_dr)
                        anchors_dr = []
                        for i in range(len(states_dr)):
                            # deriv vector will be one shorter, assume we're more interested in the end than the start
                            if i == 0:
                                anchor_dr = chspy.Anchor(times_dr[i], states_dr[i], 0.0) 
                            else:
                                anchor_dr = chspy.Anchor(times_dr[i], states_dr[i], deriv_dr[i-1])
                            anchors_dr.append(anchor_dr)
                        dr_spline = chspy.CubicHermiteSpline(n=1, anchors=anchors_dr)
                        trip_sol = dr_spline.solve(0, 
                                                   self.trip_info['limit'], 
                                                   solve_derivative = False,)
                    else:
                        solve_diff = True if self.trip_info['type'] == 'diff' else False
                        trip_sol = state.solve(self.trip_info['idx'],
                                               self.trip_info['limit'],
                                               solve_derivative = solve_diff,
                                               )
                    if trip_obj.delay:
                        self.trip_info['time'] = trip_sol[0][0] + trip_obj.delay
                    else:
                        self.trip_info['time'] = trip_sol[0][0] 
                    print(f"tripped at t = {self.trip_info['time']:.3f}")
                    print(f"state idx: {self.trip_info['idx']}")
                    print(f"limit: {tripped[1]}")
                    break
        return np.array(y)

    def equilibrium_search(self, 
                            dT, 
                            max_delay=1e10, 
                            populate_nodes=False, 
                            abs_tol=1e-10, 
                            rel_tol=1e-05, 
                            min_step = 1e-10, 
                            max_step = 10.0,
                            md_step = 1e-3,
                            abs_tol_eq = 1e-6,
                            rel_tol_eq = 1e-4,
                            max_iter = MAX_INT,
                            norm = None,
                            show_conv_metrics = False,
                            start_time = 0.0
              ):
        """
        Solves until equilibrium condition reached
        """
        self.max_delay = max_delay
        self._prepare_integrator(start_time, abs_tol, rel_tol, min_step, max_step)

        if self.trip_conditions:
            raise ValueError('equilibrium_search not compatible with trip conditions')
        
        T = []
        y = []
        y0 = np.array([self.nodes[n].y0 for n in self.nodes])
        
        diff = float('inf')
        tol = abs_tol_eq + rel_tol_eq*np.linalg.norm(y0, ord = norm)
        iters = 0
        while (diff >= tol) and (iters < max_iter):
            # find time
            if len(T) == 0:
                t_x = dT
            else:
                t_x = T[-1] + dT
            T.append(t_x)

            # calculate state 
            if (t_x <= self.max_delay):
                y.append(self.integrator.integrate_blindly(t_x, step = md_step))
            else:
                y.append(self.integrator.integrate(t_x))

            # update error & tolerance
            if len(y) == 1:
                diff = np.linalg.norm(y[-1]-y0, ord = norm)
            else:
                diff = np.linalg.norm(y[-1]-y[-2], ord = norm)
            tol = abs_tol_eq + rel_tol_eq*diff
            iters += 1

        if show_conv_metrics:
            print(f"converged at t = {T[-1]} after {iters} iterations at tol = {tol}")
            print(f"||y_k - y_{{k-1}}||_2 = {diff}")

        # populate node objects with solutions, off by default, as it can cause
        # memory blowup/leakage when running many models 
        if populate_nodes:
            self._populate_nodes(y)            

        return T, np.array(y)

    def _populate_nodes(self, sol):
        print('populating nodes objects solution vectors...')
        for s in enumerate(self.nodes.values()):
            s[1].y_out = np.array([state[s[0]] for state in sol])

    def plot_input(self, index, fac = 1.0):
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
            times, vals = np.array(times), np.array(vals)
            ax.plot(times/fac, vals)
            return ax


    def evaluate_expr(self, expr, t_array, y_array):
        """
        Evaluate a symbolic expression over time.
        
        Supports SymEngine expressions using current_y(i) by converting them to valid sympy Symbols.
        
        Parameters:
            expr: symengine.Basic or str, the symbolic expression (e.g., from p['flows'])
            t_array: (n,) array of times
            y_array: (n, m) array where y[i, j] = current_y(j) at t[i]
            
        Returns:
            (n,) array of evaluated values
        """
        # convert to sympy expression with current_y as a symbolic function
        expr_str = str(expr)
        t_sym = sp.Symbol("t")
        current_y = sp.Function("current_y")
        expr_sympy = sp.sympify(expr_str, locals={"current_y": current_y})

        # replace current_y(i) → y_i
        subs_map = {}
        indices = []

        for atom in expr_sympy.atoms(sp.Function):
            if atom.func == current_y:
                i = int(atom.args[0])
                sym = sp.Symbol(f"y_{i}")
                subs_map[atom] = sym
                indices.append(i)

        # apply the substitution
        expr_sympy = expr_sympy.subs(subs_map)

        # handle case: with or without current_y terms
        if indices:
            indices_sorted, symbols_sorted = zip(*sorted(zip(indices, subs_map.values()), key=lambda x: x[0]))
            f = sp.lambdify([t_sym] + list(symbols_sorted), expr_sympy, modules=["numpy"])
        else:
            f = sp.lambdify([t_sym], expr_sympy, modules=["numpy"])

        # evaluate expression over time
        results = []
        for i in range(len(t_array)):
            if indices:
                y_vals = y_array[i, indices_sorted]
                args = [t_array[i]] + list(y_vals)
            else:
                args = [t_array[i]]
            results.append(f(*args))

        return np.array(results)

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
        self._linked_nodes = []     # list of linked nodes
        self._dydt_linked = 0.0     # sym. expressions for linked nodes
        self.in_system = False      # flag to check if node has been added to system

    @property
    def dydt(self):
        """float: Symbolic expression for the rate of change of state variables."""
        return self.dTdt_advective + self.dTdt_internal + \
               self.dTdt_convective + self.dndt + self.dcdt + self.drdt + \
               self._dydt + self.dndt_decay + self.dydt_linked

    @dydt.setter
    def dydt(self, custom_dydt):
        # TODO: This error message shows anytime custom dydt is set, not just 
        #       when setting equal to another node. Rethink check. 
        # if isinstance(custom_dydt, Mul):
        #     print(
        #     """ Warning: You are setting this node's dynamics equal to that of 
        #         another node. If the other node's dynamics are updated, it will 
        #         not be propogated to this node. If you wish for updates to be 
        #         carried to this node, use Node.set_dydt_node() instead.
        #     """
        #         )
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
        
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
        
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
        
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

    def set_dcdt(self, 
                 n: y, beta: float, 
                 Lambda: float, 
                 lam: float, 
                 flow: bool = False, 
                 t_c: float = 0.0, 
                 t_l: float = 0.0,
                 force_steady_state: bool = False,
                 max_delay: float = 1e10):
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dcdt = 0.0
            source = n * beta / Lambda
            decay = lam * self.y()
            if flow:
                t_l = sp.Min(max_delay, t_l)
                outflow = self.y() / t_c
                if force_steady_state:
                    inflow = self.y() * sp.exp(-lam * t_l) / t_c
                else:
                    inflow = self.y(t - t_l) * sp.exp(-lam * t_l) / t_c
            else:
                inflow, outflow = 0.0, 0.0
            self.dcdt = source - decay - outflow + inflow
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
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
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
        #check that node has been added to the system
        if self.y:
            # reset in case of update
            self.dndt_decay = 0.0
            self.dndt_decay += ((n/n0) * rel_yield - lam * self.y())  
        else:
            raise ValueError("Nodes need to be added to a System() object before setting dynamics.")
        
    def set_dydt_node(self, nodes: list, coeffs: list = None):
        """
        Add dynamics of another node
        """
        if not self.in_system:
            raise ValueError('''Node dynamics cannot be set until added to 
                                a System() object''')
        if coeffs is None:
            coeffs = [1.0]*len(nodes)
        for idx, n in enumerate(nodes): 
            self._linked_nodes.append((coeffs[idx],n))

    @property
    def dydt_linked(self):
        """float: Symbolic expression for the rate of change of state variables."""
        self._dydt_linked = 0.0
        for coeff, n in self._linked_nodes:
            self._dydt_linked += (coeff*n.dydt)
        return self._dydt_linked
        

