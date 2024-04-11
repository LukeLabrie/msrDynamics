from .system_objects import Node
from jitcdde import y

class mann_HX:
    '''
    Creates nodes and sets dynamics for a fluid-to-fluid heat-exchanger using 
    a lumped (Mann's) model
    '''

    def __init__(self, m_p=0.0, c_p=0.0, w_p=0.0, p1_0=0.0, p2_0=0.0, 
                 m_s=0.0, c_s=0.0, w_s=0.0, s1_0=0.0, s2_0=0.0, m_t=0.0, 
                 c_t=0.0,t0=0.0, hA_pt=0.0, hA_st=0.0, n_lumps=1) -> None: 
        self.m_p = m_p             # primary fluid total mass
        self.c_p = c_p             # primary fluid specific heat capacity
        self.w_p = w_p             # primary fluid mass flow rate
        self.p1_0 = p1_0           # primary fluid initial inlet temp
        self.p2_0 = p2_0           # primary fluid initial outlet temp
        self.m_s = m_s             # secondary fluid total mass
        self.c_s = c_s             # secondary fluid specific heat capacity 
        self.w_s = w_s             # secondary fluid mass flow rate
        self.s1_0 = s1_0           # secondary fluid initial inlet temp
        self.s2_0 = s2_0           # secondary fluid initial inlet temp
        self.m_t = m_t             # solid total mass
        self.c_t = c_t             # solid specific heat capacity
        self.hA_pt = hA_pt         # primary->solid heat transfer coefficient 
        self.hA_st = hA_st         # secondary->solid heat transfer coefficient
        self.n_lumps = n_lumps     # number of mann models in series 
        self.nodes = None          
        self.p_inlet = None
        self.p_outlet = None
        self.s_inlet = None
        self.s_outlet = None 

        # fluid node masses
        m_pn = self.m_p/(2*self.n_lumps)
        m_sn = self.m_s/(2*self.n_lumps)

        # solid node mass 
        m_tn = self.m_t/(self.n_lumps)

        # temperatures 
        dT_p = self.p1_0-self.p2_0
        dT_s = self.s2_0-self.s1_0

        inc_p = dT_p/(2*self.n_lumps-1)
        inc_s = dT_s/(2*self.n_lumps-1)

        # define nodes
        self.nodes = {}
        for l in range(1,self.n_lumps+1):

            p1 = Node(name = f'p1_l{l}', m = m_pn, scp = self.c_p, W = self.w_p, y0 = self.p1_0-2*inc_p*l)
            p2 = Node(name = f'p2_l{l}', m = m_pn, scp = self.c_p, W = self.w_p, y0 = self.p1_0-inc_p*(l+1))
            if (l == self.n_lumps):
                self.p_outlet = p2

            t1 = Node(name = f't1_l{l}', m = m_tn, scp = self.c_t)

            s1 = Node(name = f's1_l{l}', m = m_sn, scp = self.c_s, W = self.w_s, y0 = self.s1_0+2*inc_s*l)
            s2 = Node(name = f's2_l{l}', m = m_sn, scp = self.c_s, W = self.w_s, y0 = self.s2_0+inc_s*(l+1))
            if (l == 1):
                self.s_outlet = s2

            t1.y0 = (p1.y0 * self.hA_pt + s1.y0 * self.hA_st) / (self.hA_pt + self.hA_st)

            # add nodes to instance and system 
            new_nodes = [p1,p2,t1,s1,s2]
            for n in new_nodes:
                self.nodes[n.name] = n

    def set_dynamics(self, p_inlet: y, s_inlet: y):
        '''
        Set's the dynamics of the heat exchanger nodes, per Mann's model 
        '''
        # hA per node
        hA_pt_n = self.hA_pt/(2*self.n_lumps)
        hA_st_n = self.hA_st/(2*self.n_lumps)

        for l in range(1,self.n_lumps+1):

            # get nodes of lump 
            p1 = self.nodes[f'p1_l{l}']
            p2 = self.nodes[f'p2_l{l}']
            t1 = self.nodes[f't1_l{l}']
            s1 = self.nodes[f's1_l{l}']
            s2 = self.nodes[f's2_l{l}']

            # set dynamics 
            if (l > 1):
                p1.set_dTdt_advective(source = self.nodes[f'p2_l{l-1}'].y()) 
            else:
                p1.set_dTdt_advective(source = p_inlet)
            p1.set_dTdt_convective(source = [t1.y()], hA = [hA_pt_n])

            p2.set_dTdt_advective(source = p1.y())
            p2.dTdt_convective = p1.dTdt_convective

            t1.set_dTdt_convective(source = [p1.y(),s1.y()], hA = [self.hA_pt/self.n_lumps,self.hA_st/self.n_lumps])

            if (l < self.n_lumps):
                s1.set_dTdt_advective(source = self.nodes[f's2_l{l+1}'].y())
            else:
                s1.set_dTdt_advective(source = s_inlet)
            s1.set_dTdt_convective(source = [t1.y()], hA = [hA_st_n])
            
            s2.set_dTdt_advective(source = s1.y())
            s2.dTdt_convective = s1.dTdt_convective





        
        