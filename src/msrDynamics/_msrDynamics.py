from jitcdde import jitcdde, y, t, jitcdde_input, input
import numpy as np
import chspy
import sympy as sp
import matplotlib.pyplot as plt

class System:
     '''
     Stores and operates on the system defined by a set of nodes
     '''
     def __init__(self, nodes = None, dydt = None, y0 = None) -> None:
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
          self.trip_info = None

     @property
     def dydt(self):
          return [node.dydt for node in self.nodes.values()]
     
     @dydt.setter
     def dydt(self, dydt):
          self._dydt = dydt

     def _get_full_input(self,times):
          return [f(times) for f in self.input_funcs]
     
     def _check_trip(self,time,state,diff):
          '''
          Returns true if the state variables referenced by the indicies fall outside
          of the bounds 
          '''
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
     
     
     def add_input(self, input_func):
          '''
          adds an input function of time to the integrator
          '''
          # store input function for when spline is created 
          if self.input_funcs:
               self.input_funcs.append(input_func)
          else:
               self.input_funcs = [input_func]

          # return jitcdde input object so it can be used for sympy expressions
          input_obj = input(self.n_inputs)
          self.n_inputs += 1
          return input_obj
     
     def set_custom_past(self, past: chspy.CubicHermiteSpline, t_truncate = None):
          '''
          takes a set of anchors (past) as the past for the current integrator.
          past can be truncated to t_truncate 
          '''
          if (t_truncate):
               past.truncate(t_truncate)
          self.custom_past = past

     def finalize(self, T, 
                  sdd: bool = False, 
                  md: float = 1e10,
                  max_anchors: int = 100,
                  input_tol: int = 5):
          '''
          instantializes and stores JiTCDDE integrator. Can be used for direct 
          accesss to JiTCDDE integrator object 
          T (np.ndarray): time array
          sdd (bool): state-dependent delays
          md (bool): maximum delay
          '''
          if (self.integrator): 
               self.integrator = None
          # set up system matrix
          self.dydt = [n.dydt for n in self.nodes.values()]
          if (self.input_funcs):
               # set up input spline
               print('setting up input splines...')
               spline = chspy.CubicHermiteSpline(n = self.n_inputs) 
               spline.from_function(self._get_full_input, 
                                    times_of_interest = T, 
                                    max_anchors = max_anchors, 
                                    tol = input_tol)
               self.input = spline
               # instantiate integrator
               DDE = jitcdde_input(self.dydt,self.input,max_delay = md)
               # set initial conditions
               if not (self.custom_past):
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
               DDE = jitcdde(self.dydt, max_delay = md)
               # set initial conditions
               if not (self.custom_past):
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
          '''
          Assigns state variable object y(), with unique index for each node
          '''
          index = self.n_nodes
          for n in new_nodes: 
               n.y = lambda tau=None, i = index: y(i, tau) if tau is not None else y(i)
               n.index = index
               # add node to nodes dictionary by name, if available, and by index if not
               if n.name:
                    self.nodes[n.name] = n
               else:
                    self.nodes[index] = n
               self.n_nodes += 1
               index += 1

     def get_state_by_index(self, i: int, j: int):
          '''
          returns the state of the i^th state variable at the j^th time index
          as (value, derivative)
          '''
          val = self.integrator.get_state()[j][1][i]
          deriv = self.integrator.get_state()[j][2][i]
          return (val,deriv)
               
     def solve(self, t0, tf,
               incr: float = 0.01, 
               sdd: bool = False, 
               max_delay: float = 1e10,
               max_anchors: int = 100,
               input_tol: int = 5,
               populate_nodes: bool = False,
               abs_tol: float = 1e-10,
               rel_tol: float = 1e-05):
          '''
          solves system and returns np.array() with solution matrix
          T: time array
          sdd: State-dependent delays
          max_delay: max delay
          '''
          if len(self.nodes) == 0:
               raise ValueError('No nodes have been added to the system')
          
          # clear data
          for n in self.nodes.values():
               if (n.y_out.any()):
                    n.y_out = []
          
          T = np.arange(t0, tf, incr)
          # set integrator 
          print("finalizing integrator...")
          self.finalize(T, sdd, max_delay, max_anchors, input_tol)
          self.integrator.set_integration_parameters(atol = abs_tol, rtol = rel_tol)

          # solution 
          y = []

          print("integrating...")
          if self.trip_conditions:
               self.trip_info = {}

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
                         self.trip_info['state'] = chspy.CubicHermiteSpline(n = len(self.trip_conditions), 
                                                               anchors = [chspy.Anchor(t_x, states, derivs)])

                    # check if system has tripped
                    tripped = self._check_trip(t_x,states,derivs)
                    if tripped:
                         # get trip conidition object
                         trip_obj = self.trip_conditions[tripped[0]]
                         print(y[-1][46])
                         print(f'idx {tripped[0]} tripped after integration to t = {t_x:3f} with a value of {tripped[1]}')

                         # store trip info 
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
                         trip_sol = state.solve(self.trip_info['idx'],
                                                  self.trip_info['limit'],
                                                  beginning = start,
                                                  target = self.trip_info['type'])
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
               for s in enumerate(self.nodes.values()):
                    s[1].y_out = np.array([state[s[0]] for state in y])

          return np.array(y)
     
     def plot_input(self, index):
          '''
          Returns matplotlib ax() object with plot of selected input
          '''
          if self.input is None:
               raise ValueError('No input has been defined')
          else:
               _, ax = plt.subplots()
               times = []
               vals = []
               for i in self.input:
                    times.append(i.time)
                    vals.append(i.state[index])
               ax.plot(times,vals)
               return ax

class Node:
     """
     This class represents a node in an MSR system, encapsulating physical properties and dynamics

     Attributes:
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
          y: JiTCDDE state variable object, to be assigned by the System class.
          y_out (list): Solution data, to be populated by the System class.

     Methods:
          set_dTdt_advective(source): Sets the rate of temperature change due to advection.
          set_dTdt_internal(source, k): Sets the rate of temperature change due to internal heat generation.
          set_dTdt_convective(source, hA): Sets the rate of temperature change due to convective heat transfer.
          set_dndt(r, beta_eff, Lambda, lam, C): Sets the rate of change of neutron population using the point kinetics equation.
          set_dcdt(n, beta, Lambda, lam, t_c, t_l): Sets the rate of change of precursor concentration.
          set_drdt(sources, coeffs): Sets the rate of change of reactivity based on feedback.
          dydt(): Calculates and returns the total rate of change of state variables.

     The class is designed to model the thermal-hydraulic and neutron-kinetic behavior of a node within a nuclear reactor system.
     """
     def __init__(self,
                    name: str = None,
                    m: float = 0.0,
                    scp: float = 0.0,
                    W: float = 0.0,
                    y0: float = 0.0) -> None:
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
          self._dydt = 0.0            # sym. expression for user-defined dynamics 
          self.y = None               # JiTCDDE state variable object, to be assigned by System
          self.index = None           # JiTCDDE state variable index, to be assigned by System
          self.y_out = np.array([])   # solution data, to be populated by System
          self.y_rhs = np.array([])   # solution data, to be populated by System
          self.trip_conditions = None
    
     @property
     def dydt(self):
          return self.dTdt_advective + self.dTdt_internal + \
                 self.dTdt_convective + self.dndt + self.dcdt + self.drdt + \
                 self._dydt
     
     @dydt.setter
     def dydt(self, custom_dydt):
          self._dydt = custom_dydt

     def set_dTdt_advective(self, source):
          '''
          Energy from advective heat transfer
          source: source node (node or constant)
          '''

          # check not to mix thermal and point-kinetic equations
          if self.dndt or self.dcdt or self.drdt:
               raise ValueError('''This node has already been assigned 
                                point-kinetic dynamics''')
          
          #check that node has been added to the system
          if self.y:
               # reset in case of update
               self.dTdt_advective = 0.0
               self.dTdt_advective = (source-self.y())*self.W/self.m
          else:
               raise ValueError('''Nodes need to be added to a System() object 
                                before setting dynamics.''')

     def set_dTdt_internal(self, source: callable, k: float):
          '''
          Energy from fission
          source: source of generation (state variable y(i))
          k: constant of proportionality 
          '''

          # check not to mix thermal and point-kinetic equations
          if self.dndt or self.dcdt or self.drdt:
               raise ValueError('''This node has already been assigned 
                                point-kinetic dynamics''')
          
          #check that node has been added to the system
          if self.y:
               # reset in case of update
               self.dTdt_internal = 0.0
               self.dTdt_internal = k*source/(self.m*self.scp)
          else:
               raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

     def set_dTdt_convective(self, source: list, hA: list):
          '''
          Energy from convective heat transfer
          a: (list of state variable(s) y(i)) from node(s) 
          hA: (list of state variable(s) y(i)) convective heat transfer 
               coefficient(s) * wetted area(s) (MW/C)
          '''
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
                    self.dTdt_convective += hA[i]*(source[i]-self.y())/(self.m*self.scp)
          else:
               raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

     def set_dndt(self, r: y, beta_eff: float, Lambda: float, lam: list, C: list):
          '''
          Point kinetics equation for neutron population
          r: reactivity, rho, (state variable y(i))
          beta_eff: delayed neutron fraction
          Lambda: prompt neutron lifetime (s)
          lam: decay constants for associated precursor groups
          C: precursor groups (list of state variables y(i))
          '''
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
                    precursors += g[1]*C[g[0]]
               self.dndt = (r-beta_eff)*self.y()/Lambda + precursors
          else:
               raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

     def set_dcdt(self, 
               n: y, 
               beta: float, 
               Lambda: float, 
               lam: float, 
               flow: bool = False,
               t_c: float = 0.0, 
               t_l: float = 0.0):
          '''
          Sets the derivative of the precursor concentration with respect to time (dc/dt).

          Parameters:
          - n: fractional neutron concentration (state variable y(i))
          - beta: fraction for precursor group
          - Lambda: prompt neutron lifetime
          - lam: decay constant for precursor group
          - flow: indicates if flow is considered in the model
          - t_c: core transit time (seconds)
          - t_l: loop transit time (seconds)

          Raises:
          - ValueError: if flow is True and either t_c or t_l are not set properly.
          '''
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
                    inflow = self.y(t-t_l) * sp.exp(-lam * t_l) / t_c
                    self.dcdt = source - decay - outflow + inflow
               else:
                    self.dcdt = source - decay
          else:
               raise ValueError("Nodes need to be added to a System() object before setting dynamics.")

     def set_drdt(self, sources: list, coeffs: list):
          '''
          Reactivity feedback 
          sources: list of derivatives of feedback sources (dy(i)/dt)
          coeffs: list of respective feedback coefficients
          '''
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
                    fb += s[1]*coeffs[s[0]]
               self.drdt = fb
          else:
               raise ValueError("Nodes need to be added to a System() object before setting dynamics.")