import numpy as np
import math
from collections import namedtuple

###############################################################################
# constants
###############################################################################

# domain
t0 = 0.0
tf = 1000.00
T = np.arange(t0,tf,0.01)

# reactivity insertion
inserted = 4e-3
insert_duration = 0.4/0.011 # ORNL-1845
t_ins = 300.00
t_wd = t_ins + (60*4)

# debugging constants
arbitrary_removal = 0 # test effect of increased remvoal from system 
F = 1.0               # tweak convective heat flow terms

# dimensions
pi = math.pi
R = 8.314            # ideal gas constant
m_H = 0.004          # molar mass of helium (kg/mol)
P = 2.2              # steady-state power (MW) ORNL-1845 pg. 58

# density
rho_inconel = 8.5*1000          # inconel density (kg/m^3)
rho_h = 0.167                   # helium density (kg/m^3) NEEDS TO BE TEMPERATURE DEPENDENT
rho_m = 2.75*1000               # BeO density (kg/m^3) ORNL-1845 p.

# specific heat capacities 
#scp_f = 1.9665e-3       # specific heat capacity of fuel salt (BTU/lb*defF) -> (MJ/kg-C) ORNL-TM-0728 p.8
scp_t = 0.101*4.1869e-3 # specific heat capacity of inconel 600 (BTU/lb*defF) -> (MJ/kg-C) ORNL-1845 p.113
scp_f = (0.26*4.1869e-3)  # specific heat capacity of fuel salt ORNL-1845 p.113 (BTU/lb*defF) -> (MJ/kg-C)
scp_c = 0.3*4.1869e-3   # specific heat capacity of cooolant (BTU/lb*defF) -> (MJ/kg-C) ORNL-1845 p.113
scp_h = 1.248*4.1869e-3 # specigic heat capacity of helium (BTU/lb*defF) -> (MJ/kg-C) ORNL-1845 p.113
scp_m = 0.48*4.1869e-3  # specific heat capcity of moderator (BTU/lb*defF) -> (MJ/kg-C) ORNL-1845 p.113

# delays
tau_hx_c_f = 20. # fuel-helium hx to core delay (unknown)
tau_hx_c_c = 20. # coolant-helium hx to core delay (unknown)
tau_c_hx_f = 20. # core->hx delay (unknown)
tau_h = 0.5      # helium loop delay (unknown)

# NEUTRONICS DATA
tau_l = 40.0  # ORNL-TM-0728 %16.44; % (s)
tau_c = 8.3  # ORNL-1845 p.120
#tau_l = 5.00  # ORNL-TM-0728 %16.44; % (s)
#tau_c = 8.3  # ORNL-1845 p.120
n_frac0 = 1.0  # initial fractional neutron density n/n0 (n/cm^3/s)
# Lam = 2.400E-04  # mean generation time ORNL-TM-1070 p.15 U235
Lam = (2.400E-04)  # mean generation time ORNL-TM-1070 p.15 U235
# Lam = 4.0E-04;  # mean generation time ORNL-TM-1070 p.15 U233
lam = np.array([1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00])
beta = (np.array([0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023]))  # U235
# beta = np.array([0.00023, 0.00079, 0.00067, 0.00073, 0.00013, 0.00009])  # U233
beta_t = np.sum(beta)  # total delayed neutron fraction MSRE
rho_0 = beta_t-sum(np.divide(beta,1+np.divide(1-np.exp(-lam*tau_l),lam*tau_c))) # reactivity change in going from stationary to circulating fuel
C0 = beta / Lam * (1.0 / (lam - (np.exp(-lam * tau_l) - 1.0) / tau_c))


#################################
# functions
#################################
def iwc_to_Pa(h20):
    '''
    Inches water column -> Pa
    '''
    return h20*248.8

def cfm_to_m3s(v):
    '''
    Cubic ft/min -> m3/s
    '''
    return v/2119

# fuel density (kg/m^3) ORNL-1845 pg. 113
def fuel_density(temp):
    return 1000*(4.04-0.0011*(temp+273.15))

# F to K conversion
def F_to_K(tempF):
    return (tempF-32)*5/9+273.15

def mass_flow(m,P,V,T):
    '''
    m: molar mass of gas
    P: pressure
    V: volume
    T: temperature
    '''
    return m*P*V/(R*T)       # helium mass flow (kg/s)

def h_US(C,k,D,R,r,Pr,p):
    '''
    convective heat transfer coefficient BTU/(sec*ft^2*degF) ORNL-1535 p.15
    C: leading coefficient
    k: thermal conductivity BTU/(sec*ft^2)
    D: pipe diameter (ft)
    R: Reynold's modulus
    r: Exponent for Reynold's modulus
    Pr: Prandtl modulus 
    p: Exponent for Prandtl modulus
    '''
    return 0.023*(k/D)*(R**0.8)*(Pr**0.4)

def R_US(u,r,nu):
    '''
    Reynold's Modulus ORNL-1345 p.7
    u: fluid velocity (ft/hr)
    r: pipe radius (ft)
    nu: kinematic viscosity (ft^2/hr)
    '''
    return 2*u*r/nu

def Pr_US(nu,rho,c,k):
    '''
    Prandtl Modulus ORNL-1345 p.7
    nu: kinematic viscosity (ft^2/hr)
    rho: density (lb/ft^3)
    c: heat capacity (Btu/(lb*degF))
    k: thermal conductivity BTU/(hr*ft^2)
    '''
    return nu*rho*c/k

Point = namedtuple('Point', ['x', 'y'])
def hA(W,points):
    '''
    Quadratic approximation of convective heat transfer coefficient,
    given available data

    W: flow rate (lb/s)
    points: data points of (flow rate (lb/s),hA)

    returns hA in BTU/(sec*degF)
    '''

    # Create matrices for the system of equations
    order = len(points) 
    a = [[p.x**i for i in reversed(range(order))] for p in points]
    A = np.array(a)
    b = np.array([p.y for p in points])

    # Solve for the coefficients
    coeffs = np.linalg.solve(A, b)

    # evaluate terms
    terms = [coeffs[i]*W**(order-i) for i in range(order)]

    return sum(terms)



###############################################################################
# core
###############################################################################

# wights
k_f1 = 0.475       # fractional power generation (fuel)
k_f2 = 0.475        # fractional power generation (fuel)
k_m = 1-(k_f1+k_f2) # fractional power generation (beryllium)
k_1 = 0.5
k_2 = 1-k_1

# thermal feedback (1/Kelvin, temperature provided in Kelvin) ORNL-1845 pg. 115
a_f = (-9.8e-5)*9/5
a_b = (1.1e-5)*9/5
a_c = (-5.88e-5)*9/5

# operating conditions taken from 25-hr Xenon run Exp. H-8
# temperatures 
T_fuel_avg = F_to_K(1311)       # ORNL 1845 pg. 58
T0_c_f1 = F_to_K(1209)          # core fuel inlet temp (K) ORNL-1845 pg. 120
T0_c_f2 = F_to_K(1522)          # core fuel outlet temp (K) ORNL-1845 pg. 120
T0_c_c1 = F_to_K(1226)          # core coolant inlet temp (K) ORNL-1845 pg. 121
T0_c_c2 = F_to_K(1335)          # core coolant outlet temp (K) ORNL-1845 pg. 121
T0_c_t1 = ((T0_c_f1+T0_c_c1)/2) # core tube temp
T0_c_m = F_to_K(1300)           # beryllium initial temp

# flow rates 
F_c_f = 46/15850       # core fuel flow rate (gal/min)->(m^3/s) ORNL-1845 pg. 120
F_c_c = 152/15850      # core coolant flow rate (gal/min)->(m^3/s) ORNL-1845 pg. 121

# dimensions
V_fuel = 52071.248849/1e6                                # CAD model (cm^3)->(m^3)
A_fuel = (69872.856584/2-(6*6.996))/10000                # CAD model (cm^3)->(m^2) 
V_tubes = (5453.961+34777.657)/1e6                       # CAD model (cm^3)->(m^3)
A_tubes = (73445.338-(12*7.728))/10000                   # CAD model (cm^3)->(m^2)
A_tube_bends = 2*(20044.557-(62*pi*(3.137/2)**2))/10000  # CAD model (cm^2)->(m^2) 
V_coolant = (248933.207)/1e6                             # CAD model (cm^3)->(m^3)
A_mc = 961677.131/10000                                  # CAD model (cm^2)->(m^2)
V_m = 926899.473/1e6                                     # CAD model (cm^2)->(m^3)

# density 
# rho_c = 1000*0.78  # coolant density (kg/m^3) 

# density data from https://apps.dtic.mil/sti/tr/pdf/AD0622191.pdf
temperature = np.array([F_to_K(t) for t in [1576.8, 1886.2, 2093.5, 2268.9, 2491.2]])
density = np.array([16.0186*r for r in [46.738, 43.926, 42.224, 40.641, 38.933]]) # lb/ft^3 -> kg/m^3

# line of best fit 
coefficients = np.polyfit(temperature, density, 1)
slope, intercept = coefficients
def coolant_density(T_c):
    return T_c*slope + intercept
rho_c = coolant_density((T0_c_c1+T0_c_c2)/2)

# mass
m_f_c = fuel_density(T_fuel_avg)*V_fuel
m_c_c = rho_c*V_coolant         # coolant mass (kg)
m_m_c = (5490/2.205)            # ORNL-1845 p.111s (lb)->(kg)
m_t     = V_tubes*rho_inconel   # mass of tubes (kg)  

# mass flow rate 
W_f = F_c_f * fuel_density(T_fuel_avg)   # fuel mass flow rate (kg/s)
W_c = F_c_c * rho_c                      # coolant mass flow rate (kg/s)

# convective heat transfer coefficient for the core fuel
rho_c_f = 187                           # core fuel density in (lb/ft^3) ORNL-1535 p.16
cp_c_f = 0.26                           # core fuel heat capacity (Btu/(lb*degF)) ORNL-1535 p.16
mu_c_f = 60*60*8.27e-3                  # absolute viscosity of core fuel lb/(hr*ft) ORNL-1535 p.16 
nu_c_f = mu_c_f/rho_c_f                 # kinematic viscosity (ft^2/hr)
d_c_f = 0.097933071                     # fuel pipe diameter (ft) CAD Model
v_c_f = (40*8.02083)/(pi*(d_c_f/2)**2)  # velocity of core fuel (ft/hr) (gpm->ft^3/hr conversion)
k_f_US = 4.17e-4                        # thermal conductivity of the fuel in BTU/(sec*ft^2) ORNL-1535 p.16
k_f_US_hr = 60*60*4.17e-4               # thermal conductivity of the fuel in BTU/(hr*ft^2) ORNL-1535 p.16
R_c_f = R_US(v_c_f,d_c_f/2,nu_c_f)      # Reynold's number for core fuel 
r_c_f = 0.8                             # coefficient on Reynold's modulus
Pr_c_f = Pr_US(nu_c_f,rho_c_f,cp_c_f,k_f_US_hr) # Prandtl modulus core fuel
p_c_f = 0.4                                     # Coefficient on Prandtl's modulus
C_c_f = 0.023                                   # leading coefficient
h_f_US = h_US(C_c_f,k_f_US,d_c_f,R_c_f,r_c_f,Pr_c_f,p_c_f) # Convective heat transfer coefficient core fuel (BTU/(sec*ft^2*degF))
h_f_US_hr = h_f_US*60*60                        # Convective heat transfer coefficient core fuel (BTU/(hr*ft^2*degF))
h_f_c = ((h_f_US_hr/0.1761)*1e-6)                   # Convective heat transfer coefficient core fuel (MW/(m^2*degC))
hA_f_c = h_f_c*A_fuel

# convective heat transfer coefficient for the core tubes
hA_t_hx_US = 1/0.137 # BTU/(sec*degF) ORNL-1535 p. 47
h_t_US = (hA_t_hx_US/11.15+hA_t_hx_US/8.02)/2 # tube htc (Btu/(sec*ft^2*defF))
hA_t_c_US = h_t_US*A_fuel*10.764
hA_t_c = hA_t_c_US*(9/5)*(1.05504)*(1e-3)  # MW/(degK)

# fuel-tube coefficient (core)
hA_ft_c = 1/((1/hA_t_c)+(1/hA_f_c))

# coolant
A_coolant_tubes_US = 317085/929      # CAD model (cm^2)->(ft^2)
h_c_c_US = 166                       # coolant heat transfer coefficient (BTU/(hr*ft^2*defF)) ORNL-1535 p.23
hA_c_c_US = A_coolant_tubes_US*h_c_c_US

# elbows
hA_t_c12_US = 1/(0.130)                       # ORNL-1535 pg.24
hA_tc_c_US = 1/((1/hA_c_c_US)+(1/hA_t_c12_US)) 
hA_tc_c = hA_tc_c_US*(9/5)*(1.05504)*(1e-3)   # MW/(degK)

# moderator
hA_m_US = 1/2.060 # ORNL-1535 p.28
hA_c_US = 1000*1/0.771 # ORNL-1535 p.28
hA_mc_US = 1/((1/hA_m_US)+(1/hA_c_US))
hA_mc_c = hA_mc_US*(9/5)*(1.05504)*(1e-3)  # MW/(degK)

# specific heat 
mcp_t_c = scp_t*m_t
mcp_f_c = scp_f*m_f_c
mcp_c_c = scp_c*m_c_c 
mcp_m_c = scp_m*m_m_c

###############################################################################
# fuel-helium heat exchanger
###############################################################################

# initial temperatures
T0_hfh_f1 = F_to_K(1450)   # fuel-helium hx fuel inlet temp (K) ORNL-1845 pg. 121
T0_hfh_f2 = F_to_K(1150)   # fuel-helium hx fuel oulet temp (K) ORNL-1845 pg. 121
T0_hfh_h1 = F_to_K(180)    # fuel-helium hx helium inlet temp (K) ORNL-1845 pg. 121
T0_hfh_h2 = F_to_K(620)    # fuel-helium hx helium outlet temp (K) ORNL-1845 pg. 121
T0_hfh_t1 = ((T0_hfh_f1+T0_hfh_h1)/2) # initial tube temp

# flow rates 
F_hfh_h1 = (7300)/2119    # fuel-helium hx helium fuel inlet flow rate (ft^3/min)->(m^3/s) ORNL-1845 pg. 121
F_hfh_h2 = (12300)/2119    # fuel-helium hx helium fuel outlet flow rate (ft^3/min)->(m^3/s) ORNL-1845 pg. 121

# dimensions
L_eff_US = 93.65*12  # effective tube length of hx (in) ORNL-1535 pg. 47
V_p_hx_US = (math.pi*((1.0/2)-0.109)**2)*L_eff_US # volume in hx tubes (in^3) ORNL-1535 pg. 47
V_p_hx = V_p_hx_US/61020 # in^3 -> m^3
V_t_hx_US = ((L_eff_US)*math.pi*((1.0/2))**2 - V_p_hx_US) # tube volume in heat exchangers (in^3)
V_t_hx = V_t_hx_US/61020 # in^3 -> m^3
A_t_hx = (pi*((1.0)-2*0.109)*L_eff_US)/144 # hx inner tube area (in^2 -> ft^2) ORNL-1535 p.47  
A_to_hx = (pi*(1.0/2)*93.65*12)/144 # hx outer tube area (in^2 -> ft^2) ORNL-1535 p.47  

# mass 
m_f_hx = V_p_hx*fuel_density(T_fuel_avg)
m_h_hxfh = (((27.0*27.5*27)/61020)-V_t_hx-V_p_hx)*rho_h # hx volume minus tube volume
m_t_hxfh = V_t_hx*rho_inconel

# heat transfer fuel<->tube
hA_f_hx_US = 5.724   # BTU/(sec*degF) ORNL-1535 p.47
hA_t_hx_US = 1/0.137 # BTU/(sec*degF) ORNL-1535 p.47
hA_ft_hx_US = 1/((1/hA_f_hx_US)+(1/hA_t_hx_US))
hA_ft_hx = hA_ft_hx_US*(9/5)*(1.05504)*(1e-3) 

# mass flow rate of helium
W_h_fh = F_hfh_h1*(0.1362) # (kg/s) = volumetric flow rate * denisty of helium at 180F, 2.2 H20, 7300 cfm ORNL-1845 p.122
W_h_fh_US = W_h_fh*2.205   # (lb/s)

# available coefficient data for helium ORNL-1535 p.47 
fh_p1 = Point(0.5,0.58)
fh_p2 = Point(1.0,0.8)
fh_p3 = Point(1.5,0.955)
hA_h_hx_US = hA(W_h_fh_US,[fh_p1,fh_p2,fh_p3])
hA_ht_hx_US = 1/((1/hA_h_hx_US)+(1/hA_t_hx_US))
hA_ht_hx = hA_ht_hx_US*(9/5)*(1.05504)*(1e-3) # BTU/(sec*degF) -> (MW/C)

# product of mass and specific heat capacities 
mcp_t_hx = m_t_hxfh*scp_t
mcp_f_hx = scp_f*m_f_hx
mcp_h_hxfh = m_h_hxfh*scp_h

###############################################################################
# helium-water heat exchanger (fuel loop)
###############################################################################

# initial temperatures
T0_hhwf_h1 = F_to_K(620)    # helium-water hx (fuel loop) helium inlet temp (K) ORNL-1845 pg. 121
T0_hhwf_h2 = F_to_K(180)    # helium-water hx (fuel loop) helium outlet temp (K) ORNL-1845 pg. 121

T0_hhwf_w1 = F_to_K(61)    # helium-water hx (fuel loop) water inlet temp (K) ORNL-1845 pg. 121 
T0_hhwf_w2 = F_to_K(124)   # helium-water hx (fuel loop) water water outlet temp (K) ORNL-1845 pg. 121

T0_hhwf_t1 = ((T0_hhwf_h1+T0_hhwf_w1)/2)

V_h_hxhw_US = (pi*((0.625/2)-0.049)**2)*825*12 # in^3 ORNL-1535 p.47
V_h_hxhw = V_h_hxhw_US/61020 # m^3 
m_h_hxhwf = V_h_hxhw*rho_h

V_t_hxhw_US = ((pi*((0.625/2))**2)*825*12)-V_h_hxhw_US
V_t_hxhw = V_t_hxhw_US/61020 # m^3
m_t_hxhwf = rho_inconel*V_t_hxhw
mcp_t_hxhw = m_t_hxhwf*scp_t

# water mass flow rate
W_hhwf_w = 998*((103*2)/15850)                 # water flow (kg/s) ORNL-1845 p.121
W_hhwf_w_US = W_hhwf_w*2.205                 # (lb/s)
V_w_US = (27*27.5*27)-V_t_hx_US-V_h_hxhw_US  # in^3
V_w = V_w_US/61020                           # in^3 -> m^3
m_w_hxhwf = V_w*998 
scp_w = 4.181e-3
mcp_w = m_w_hxhwf*scp_w
mcp_h_hxhw = m_h_hxhwf*scp_h

# helium<->tube
hxhw_h_p1 = Point(0.5,1.40)
hxhw_h_p2 = Point(1.0,2.37)
hxhw_h_p3 = Point(1.5,3.22)
hA_h_hxhw_US = hA(W_h_fh_US,[hxhw_h_p1,hxhw_h_p2,hxhw_h_p3])
hA_t_hxhw_US = 1/0.00456 # BTU/(sec*degF)             # ORNL-1535 p.47
hA_ht_hxhw_US = 1/((1/hA_h_hxhw_US)+(1/hA_t_hxhw_US))
hA_ht_hxhw = hA_ht_hxhw_US*(9/5)*(1.05504)*(1e-3)     # BTU/(sec*degF) -> MW/C

# tube<->water
hxhw_w_p1 = Point(3.54,12.5)
hxhw_w_p2 = Point(7.08,22.8)
hxhw_w_p3 = Point(10.61,30.8)
hA_w_hxhw_US = hA(W_hhwf_w_US,[hxhw_w_p1,hxhw_w_p2,hxhw_w_p3])
hA_tw_hxhw_US = 1/((1/hA_t_hxhw_US)+(1/hA_w_hxhw_US))
hA_tw_hxhw = hA_tw_hxhw_US*(9/5)*(1.05504)*(1e-3)     # BTU/(sec*degF) -> MW/C

###############################################################################
# coolant-helium heat exchanger 
###############################################################################

# initial temperatures
T0_hch_c1 = F_to_K(1235)  # coolant-helium hx coolant inlet temp (K) ORNL-1845 pg. 121 
T0_hch_c2 = F_to_K(1105)  # coolant-helium hx coolant outlet temp (K) ORNL-1845 pg. 121 

T0_hch_h1 = F_to_K(170)  # coolant-helium hx helium inlet temp (K) ORNL-1845 pg. 122 
T0_hch_h2 = F_to_K(1020)  # coolant-helium hx helium outlet temp (K) ORNL-1845 pg. 122 

T0_hch_t1 = ((T0_hch_c1+T0_hch_h1)/2)

# dimensions
V_c_hxch_US = 101.6*12*pi*((1.0/2)-0.109)**2 # (in^3) ORNL-1535 p.58
V_c_hxch = V_c_hxch_US/61020 # (in^3) -> (m^3)
V_t_hxch_US = (101.6*12*pi*((1.0/2))**2)-V_c_hxch_US # (in^3) ORNL-1535 p.58
V_t_hxch = V_t_hxch_US/61020 # (in^3) -> (m^3)
V_h_hxch_US = (16.25**2)*13.75 # (in^3) ORNL-1535 p.58
V_h_hxch = V_h_hxch_US/61020 # (in^3) -> (m^3)

# mass
m_c_hx = V_c_hxch * rho_c
m_t_hxch = V_t_hxch*rho_inconel
m_h_hxch = rho_h*V_h_hxch 

# product of mass and specific heat capacity
mcp_c_hxch = m_c_hx*scp_c    
mcp_h_hxch = m_h_hxch*scp_h
mcp_t_hxch = m_t_hxch*scp_t

# helium mass flow rate 
F_h_ch = 2000/2119 # ft^3/min->m^3/s ORNL-1845 p.122
W_h_ch = F_h_ch*0.1389 # volumetric flow rate * denisty of helium at 170F, 4.0 H20 ORNL-1845 p.122
W_hxch_h_US = W_h_ch*2.205                 # (kg/s)->(lb/s)

# coolant<->tube
hA_c_hx_US = 17.3 # (BTU/(sec*degF)) ORNL-1535 p.58
hA_t_hxch_US = 1/0.126 # (BTU/(sec*degF)) ORNL-1535 p.58
hA_ct_US_hx = 1/((1/hA_c_hx_US)+(1/hA_t_hxch_US))  # BTU/(sec*degF)
hA_ct_hx = hA_ct_US_hx*(9/5)*(1.05504)*(1e-3) # MW/C

# tube<->helium
hxch_h_p1 = Point(0.11,0.526)
hxch_h_p2 = Point(0.22,0.743)
hxch_h_p3 = Point(0.33,0.889)
hA_h_hxch_US = hA(W_hxch_h_US,[hxch_h_p1,hxch_h_p2,hxch_h_p3])  
hA_th_hxch_US = 1/((1/hA_h_hxch_US)+(1/hA_t_hxch_US))
hA_th_hxch = hA_th_hxch_US*(9/5)*(1.05504)*(1e-3) # BTU/(sec*degF) -> MW/C

###############################################################################
# helium-water heat exchanger (coolant loop) 
###############################################################################

# initial temperatures 
T0_hhwc_h1 = F_to_K(1020)    # helium-water hx (coolant loop) helium inlet temp (K) ORNL-1845 pg. 122
T0_hhwc_h2 = F_to_K(170)    # helium-water hx (coolant loop) helium outlet temp (K) ORNL-1845 pg. 122

T0_hhwc_w1 = F_to_K(70)    # helium-water hx (coolant loop) water inlet temp (K) ORNL-1845 pg. 121 
T0_hhwc_w2 = F_to_K(100)   # helium-water hx (coolant loop) water water outlet temp (K) ORNL-1845 pg. 121

T0_hhwc_t1 = ((T0_hhwc_h1+T0_hhwc_w1)/2)

# dimensions
V_h_hxhwc_US = (pi*((0.625/2)-0.049)**2)*255.1*12 # in^3 ORNL-1535 p.58
V_h_hxhwc = V_h_hxhwc_US/61020 # (in^3) -> (m^3)
V_t_hxhwc_US = ((pi*((0.625/2))**2)*255.1*12)-V_h_hxhwc_US # in^3 ORNL-1535 p.58
V_t_hxhwc = V_t_hxhwc_US/61020 # (in^3) -> (m^3)
V_w_hxhwc_US = (17**3)-V_t_hxhwc_US-V_h_hxhwc_US # (in^3) ORNL-1535 p.58
V_w_hxhwc = V_w_hxhwc_US/61020 # (in^3) -> (m^3)

# mass
m_h_hxhwc = rho_h*V_h_hxhwc
m_t_hxhwc = V_t_hxhwc*rho_inconel
m_w_hxhwc = V_w_hxhwc*998

# product of mass and specific heat capacity
mcp_h_hxhwc = m_h_hxhwc*scp_h
mcp_t_hxhwc = m_t_hxhwc*scp_t
mcp_w_hxhwc = m_w_hxhwc*scp_w

# tube<->helium
hxhwc_h_p1 = Point(0.11,0.424)
hxhwc_h_p2 = Point(0.22,0.709)
hxhwc_h_p3 = Point(0.33,0.957)
hA_h_hxhwc_US = hA(W_hxch_h_US,[hxhwc_h_p1,hxhwc_h_p2,hxhwc_h_p3])
hA_t_hxhwc_US = 1/0.015 # ORNL-1535 p.58
hA_ht_hxhwc_US = 1/((1/hA_h_hxhwc_US)+(1/hA_t_hxhwc_US))
hA_ht_hxhwc =  hA_ht_hxhwc_US*(9/5)*(1.05504)*(1e-3) # BTU/(sec*degF) -> MW/C

# mass flow 
W_hhwc_w = 998*((38.3*2)/15850)  # water flow rate in helium-water hx (kg/s) ORNL-1845 p.122
W_hhwc_w_US = W_hhwc_w*2.205 # (kg/s)->(lb/s)

# tube<->water
hxhwc_w_p1 = Point(3.33,5.92)
hxhwc_w_p2 = Point(6.67,10.2)
hxhwc_w_p3 = Point(10.0,14.3)
hxhwc_w_p4 = Point(16.67,21.5)
hA_w_hxhwc_US = hA(W_hhwc_w_US,[hxhwc_w_p1,hxhwc_w_p2,hxhwc_w_p3,hxhwc_w_p4])
hA_tw_hxhwc_US = 1/((1/hA_t_hxhwc_US)+(1/hA_w_hxhwc_US))
hA_tw_hxhwc = hA_tw_hxhwc_US*(9/5)*(1.05504)*(1e-3) # BTU/(sec*degF) -> MW/C
















