from jitcdde import y, t
import numpy as np

class Node:
    def __init__(self,
                 m: float = 0.0,
                 scp: float = 0.0,
                 W: float = 0.0,
                 y0: float = 0.0) -> None:
        self.m = m              # mass (kg)
        self.scp = scp          # specific heat capacity (J/(kg*K))
        self.W = W              # mass flow rate (kg/s)
        self.y0 = y0            # initial temperature (K)
        self.BulkFlow = 0.0
        self.internal = 0.0
        self.convective = 0.0
        self.dndt = 0.0
        self.dcdt = 0.0
        self.drdt = 0.0
        self.y = None           # temperature (K), to be assigned by System class
    
    def set_dTdt_bulkFlow(self, source):
        '''
        Energy from bulk flow
        source: source node (state variable y(i) or constant)
        dumped: if 'from node' is a constant (indicates dumping instead of 
                recirculation), this needs to be set to true
        '''
        self.BulkFlow = (source-self.y())*self.W/self.m
 

    def set_dTdt_internal(self, source: callable, k: float):
        '''
        Energy from fission
        source: source of generation (state variable y(i))
        k: constant of proportionality 
        '''
        self.internal = k*source
    

    def set_dTdt_convective(self, source: list, hA: list):
        '''
        Energy from convective heat transfer
        a: from node(s) (state variable(s) y(i))
        hA_mcp: ratio of [convective heat transfer coefficient(s) * wetted area(s) (MW/C)]
        '''
        for i in range(len(source)):
                self.convective += hA[i]*(source[i]-self.y())/(self.m*self.scp)
    
    def set_dndt(self, r: y, beta_eff: float, Lambda: float, lam: list, C: list):
         '''
         Point kinetics equation for neutron population
         r: reactivity, rho, (state variable y(i))
         beta_eff: delayed neutron fraction
         Lambda: prompt neutron lifetime (s)
         lam: decay constants for associated precursor groups
         C: precursor groups (list of state variables y(i))
         '''
         precursors = 0.0
         for g in enumerate(lam):
              precursors += g[1]*C[g[0]]
         self.dndt = (r-beta_eff)*self.y()/Lambda + precursors

    def set_dcdt(self, 
                 n: y, 
                 beta: float, 
                 Lambda: float, 
                 lam: float, 
                 t_c: float, 
                 t_l: float):
         '''
         Precursor concentration
         n: fractional neutron concentration (state variable y(i))
         beta: fraction for precuror group
         Lambda: prompt neutron lifetime 
         lam: decay constant for precursor group
         t_c: core transit time (s)
         t_l: loop transit time (s)
         '''
         source = n*beta/Lambda
         decay = lam*self.y()
         outflow = self.y()/t_c
         inflow = self.y(t-t_l)*np.exp(-lam*t_l)/t_c
         self.dcdt = source - decay - outflow + inflow 
    
    def set_drdt(self, sources: list, coeffs: list):
         '''
         Reactivity feedback 
         sources: list of derivatives of feedback sources (dy(i)/dt)
         coeffs: list of respective feedback coefficients
         '''
         fb = 0.0
         for s in enumerate(sources):
              fb += s[1]*coeffs[s[0]]
         self.drdt = fb
    
    def dydt(self):
          y1 = self.BulkFlow
          y2 = self.internal
          y3 = self.convective
          y4 = self.dndt
          y5 = self.dcdt
          y6 = self.drdt
          sum = y1 + y2 + y3 + y4 + y5 + y6
          return sum

class System:
     def __init__(self) -> None:
          self.nNodes = 0
          self.nodes = []

     def addNodes(self, newNodes: list):
          index = self.nNodes
          for n in newNodes: 
                n.y = lambda tau=None, i = index: y(i, tau) if tau is not None else y(i)
                self.nodes.append(n)
                self.nNodes += 1
                index += 1
          return None