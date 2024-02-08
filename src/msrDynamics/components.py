from msrDynamics.objects import Node

class mann_HX:
    '''
    Creates nodes and sets dynamics for a fluid-to-fluid heat-exchanger using 
    a lumped (Mann's) model
    '''

    def __init__(self, system, m_p=0.0, c_p=0.0, w_p=0.0, p1_0=0.0, p2_0=0.0, m_s=0.0, 
                 c_s=0.0, w_s=0.0, s1_0=0.0, s2_0=0.0, m_t=0.0, c_t=0.0,
                 t0=0.0, hA_pt=0.0, hA_st=0.0, n_lumps=1, p_inlet = None, 
                 s_inlet = None) -> None:
        self.system = system       # msrDynamics System object
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
        self.nodes = None          # Nodes, populated by get_nodes

    def get_nodes(self):

        nodes = []

        # fluid node masses
        m_pn = self.m_p/(2*self.n_lumps)
        m_sn = self.m_s/(2*self.n_lumps)

        # solid node mass 
        m_tn = self.m_t/(self.n_lumps)

        # hA per node
        hA_pt_n = self.hA_pt/(2*self.n_lumps)
        hA_st_n = self.hA_st/(2*self.n_lumps)

        # temperatures 
        dT_p = self.p1_0-self.p2_0
        dT_s = self.s2_0-self.s1_0

        inc_p = dT_p/(2*self.n_lumps-1)
        inc_s = dT_s/(2*self.n_lumps-1)

        # instantiate nodes
        for l in range(self.n_lumps):

            p1 = Node(m = m_pn, scp = self.c_p, W = self.w_p, y0 = self.p1_0-2*inc_p*l)
            p2 = Node(m = m_pn, scp = self.c_p, W = self.w_p, y0 = self.p1_0-inc_p*(l+1))

            t1 = Node(m = m_tn, scp = self.c_t)

            s1 = Node(m = m_sn, scp = self.c_s, W = self.w_s, y0 = self.s1_0+2*inc_s*l)
            s2 = Node(m = m_sn, scp = self.c_s, W = self.w_s, y0 = self.s2_0+inc_s*(l+1))

            t1.y0 = (p1.y0 * hA_pt_n + s1.y0 * hA_st_n) / (hA_pt_n + hA_st_n)

            # add nodes to system 
            new_nodes = [p1,p2,t1,s1,s2]
            self.system.add_nodes(new_nodes)

            for n in new_nodes:
                nodes.append(n)

            # store nodes in instance
            if self.nodes is None:
                self.nodes = new_nodes
            else:
                for n in new_nodes:
                    self.nodes.append(n)

        # define dynamics 
        for l in range(self.n_lumps):
            p1 = self.nodes[0+l*5]
            p2 = self.nodes[1+l*5]
            t1_p = self.nodes[2+l*5]
            t1_s = self.nodes[2+(self.n_lumps-l-1)*5]
            s1 = self.nodes[3+l*5]
            s2 = self.nodes[4+l*5]
            if (self.n_lumps > 1) and (l > 0):
                p1.set_dTdt_advective(source = nodes[-(9*l)].y()) 
            p1.set_dTdt_convective(source = [t1_p.y()], hA = [hA_pt_n])

            p2.set_dTdt_advective(source = p1.y())
            p2.dTdt_convective = p1.dTdt_convective

            t1_p.set_dTdt_convective(source = [p1.y(),nodes[-(5*l+2)].y()], hA = [2*hA_pt_n,2*hA_st_n])
            # t1_s.set_dTdt_convective(source = [p1.y(),s1.y()], hA = [2*hA_pt_n,2*hA_st_n])

            if (self.n_lumps > 1) and (l > 0):
                s1.set_dTdt_advective(source = nodes[-(5*l+1)].y())
            s1.set_dTdt_convective(source = [t1_s.y()], hA = [hA_st_n])
            
            s2.set_dTdt_advective(source = s1.y())
            s2.dTdt_convective = s1.dTdt_convective


        self.nodes = nodes
        return self.nodes





        
        