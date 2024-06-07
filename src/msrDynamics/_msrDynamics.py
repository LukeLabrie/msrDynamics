from jitcdde import jitcdde, y, t, jitcdde_input, input
import numpy as np
from chspy import CubicHermiteSpline
from chspy._chspy import solve_from_anchors
import sympy as sp

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
          self.trip_conditions = None
          self.custom_past = None

     @property
     def dydt(self):
          return [node.dydt for node in self.nodes.values()]
     
     def _get_full_input(self,times):
          return [f(times) for f in self.input_funcs]
     
     def _check_trip(self,state,diff,bounds,trip_types):
          '''
          Returns true if the state variables referenced by the indicies fall outside
          of the bounds 
          '''
          # print(state)
          # print(diff)
          # print(bounds)
          # print(trip_types)
          for idx in range(len(state)):
               bounds_s = bounds[idx]
               min_s = bounds_s[0]
               max_s = bounds_s[1]
               if trip_types[idx] == 'state':
                    s = state[idx]
               elif trip_types[idx] == 'diff':
                    try:
                         s = diff[idx]
                    except:
                         return None
               else:
                    raise ValueError('''Invalid trip type. Currently supported 
                                        are 'state' and 'diff'.''')
               # print(s)
               if (s < min_s):
                    return (idx, min_s)
               elif (s > max_s):
                    return (idx, max_s)
          return None
     
     @dydt.setter
     def dydt(self, dydt):
          self._dydt = dydt
     
     def add_input(self, input_func, T):
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
     
     def set_trip_conditions(self, indices, bounds, trip_types = None):
          '''
          set bounds for state variables given by index in indicies.
          when bound is reached, integration stops
          '''
          self.trip_conditions = {}
          self.trip_conditions['indicies'] = indices
          self.trip_conditions['bounds'] = bounds 
          if trip_types:
               self.trip_conditions['trip_types'] = trip_types
          else:
               self.trip_conditions['trip_types'] = ['state']*len(indices)
     
     def set_custom_past(self, past: CubicHermiteSpline, t_truncate = None):
          '''
          takes a set of anchors (past) as the past for the current integrator.
          past can be truncated to t_truncate 
          '''
          if (t_truncate):
               past.truncate(t_truncate)
          self.custom_past = past

     def finalize(self, T, 
                  sdd: bool = False, 
                  md: float = 1e10):
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
               spline = CubicHermiteSpline(n = self.n_inputs) 
               spline.from_function(self._get_full_input, times_of_interest = T)
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
               
     def solve(self, 
               T: list, 
               sdd: bool = False, 
               max_delay: float = 1e10):
          '''
          solves system and returns np.array() with solution matrix
          T: time array
          sdd: State-dependent delays
          max_delay: max delay
          '''
          # clear data
          for n in self.nodes.values():
               if (n.y_out.any()):
                    n.y_out = []

          # set integrator 
          print("finalizing integrator...")
          self.finalize(T, sdd, max_delay)

          # solution 
          y = []

          print("integrating...")
          if self.trip_conditions:
          # integrate with trip conditions
               for t_x in T:
                    y.append(np.array(self.integrator.integrate(t_x)))
                    states = y[-1][self.trip_conditions['indicies']]
                    try:
                         derivs = (y[-1]-y[-2])/t_x
                    except:
                         derivs = []
                    tripped = self._check_trip(states,
                                               derivs,
                                               self.trip_conditions['bounds'],
                                               self.trip_conditions['trip_types'])
                    if tripped:
                         self.trip_conditions['tripped_idx'] = self.trip_conditions['indicies'][tripped[0]]
                         self.trip_conditions['tripped_val'] = tripped[1]
                         print(f"tripped at t = {t_x}: ({tripped[0]},{tripped[1]})")
                         # calculate trip time
                         window = 1
                         trip_sol = []
                         while not trip_sol:
                              sol_temp = solve_from_anchors((self.integrator.get_state()[-1-window],
                                                            self.integrator.get_state()[-1]),
                                                            self.trip_conditions['tripped_idx'],
                                                            self.trip_conditions['tripped_val'])
                              trip_sol = sol_temp
                              window += 1
                         self.trip_conditions['t_trip'] = trip_sol[0][0]
                         break
          else:
               for t_x in T:
                    y.append(self.integrator.integrate(t_x))

          # populate node objects with solutions 
          for s in enumerate(self.nodes.values()):
               s[1].y_out = np.array([state[s[0]] for state in y])

          return np.array(y)

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