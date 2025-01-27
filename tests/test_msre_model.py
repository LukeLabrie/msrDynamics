from _msre_parameters import get_params
import numpy as np
from jitcdde import t
from msrDynamics import Node, System
import os
import pandas as pd
import numpy.typing as npt

# msre model builder using msrDynamics
def _build_msre_model(power: float, times: npt.ArrayLike, t_ins: float = 2500, pcm:float = 5e-4):

    p = get_params(power)

    MSRE = System()

    # radiator
    T_out_rc = Node(m=p["mn_rp"], scp=p["mcp_rpn"] / p["mn_rp"], W=p["W_rp"], y0=p["T0_rp"], name = 'T_out_rc')
    T_out_air = Node(m=p["mn_rs"], scp=p["mcp_rsn"] / p["mn_rs"], W=p["W_rs"], y0=p["T0_rs"], name='T_out_air')

    # heat exchanger
    T_hf1 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p1"], name='T_hf1')
    T_hf2 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p2"], name='T_hf2')
    T_hf3 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p3"], name='T_hf3')
    T_hf4 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p4"], name='T_hf4')
    T_ht1 = Node(m=p["m_tn"], scp=p["scp_t"], y0=p["T0_t1"], name='T_ht1')
    T_ht2 = Node(m=p["m_tn"], scp=p["scp_t"], y0=p["T0_t2"], name='T_ht2')
    T_hc1 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s1"], name='T_hc1')
    T_hc2 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s2"], name='T_hc2')
    T_hc3 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s3"], name='T_hc3')
    T_hc4 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s4"], name='T_hc4')

    # core
    n = Node(y0=p["n_frac0"], name='n')
    C1 = Node(y0=p["C0"][0], name='C1')
    C2 = Node(y0=p["C0"][1], name='C2')
    C3 = Node(y0=p["C0"][2], name='C3')
    C4 = Node(y0=p["C0"][3], name='C4')
    C5 = Node(y0=p["C0"][4], name='C5')
    C6 = Node(y0=p["C0"][5], name='C6')
    rho = Node(y0=0.0, name='rho')

    # add reactivity input
    def rho_insert(t):
        if t < t_ins:
            return 0.0
        else:
            return pcm

    rho_ext = MSRE.add_input(rho_insert, times)

    T_cg = Node(m=p["mcp_g1"] / p["scp_g"], scp=p["scp_g"], y0=p["T0_g1"], name='T_cg')
    T_cf1 = Node(m=p["mn_f"], scp=p["scp_f"], W=p["W_f"], y0=p["T0_f1"], name='T_cf1')
    T_cf2 = Node(m=p["mn_f"], scp=p["scp_f"], W=p["W_f"], y0=p["T0_f2"], name='T_cf2')

    MSRE.add_nodes([
        T_out_rc, T_out_air, T_hf1, T_hf2, T_hf3, T_hf4, T_ht1, T_ht2,
        T_hc1, T_hc2, T_hc3, T_hc4, n, C1, C2, C3, C4, C5, C6, T_cg, T_cf1, T_cf2, rho
    ])

    # dynamics
    # radiator
    T_out_rc.set_dTdt_advective(source=T_hc4.y(t - p["tau_hx_r"]))
    T_out_rc.set_dTdt_convective(source=[T_out_air.y()], hA=[p["hA_rpn"]])

    T_out_air.set_dTdt_advective(source=p["Trs_in"])
    T_out_air.set_dTdt_convective(source=[T_out_rc.y()], hA=[p["hA_rsn"]])

    # heat exchanger
    T_hf1.set_dTdt_advective(source=T_cf2.y(t - p["tau_c_hx"]))
    T_hf1.set_dTdt_convective(source=[T_ht1.y()], hA=[p["hA_pn"]])

    T_hf2.set_dTdt_advective(source=T_hf1.y())
    T_hf2.dTdt_convective = T_hf1.dTdt_convective

    T_hf3.set_dTdt_advective(source=T_hf2.y())
    T_hf3.set_dTdt_convective(source=[T_ht2.y()], hA=[p["hA_pn"]])

    T_hf4.set_dTdt_advective(source=T_hf3.y())
    T_hf4.dTdt_convective = T_hf3.dTdt_convective

    T_ht1.set_dTdt_convective(
        source=[T_hf1.y(), T_hf1.y(), T_hc3.y(), T_hc3.y()],
        hA=[p["hA_pn"], p["hA_pn"], p["hA_sn"], p["hA_sn"]]
    )
    T_ht2.set_dTdt_convective(
        source=[T_hf3.y(), T_hf3.y(), T_hc1.y(), T_hc1.y()],
        hA=[p["hA_pn"], p["hA_pn"], p["hA_sn"], p["hA_sn"]]
    )

    T_hc1.set_dTdt_advective(source=T_out_rc.y(t - p["tau_r_hx"]))
    T_hc1.set_dTdt_convective(source=[T_ht2.y()], hA=[p["hA_sn"]])

    T_hc2.set_dTdt_advective(source=T_hc1.y())
    T_hc2.dTdt_convective = T_hc1.dTdt_convective

    T_hc3.set_dTdt_advective(source=T_hc2.y())
    T_hc3.set_dTdt_convective(source=[T_ht1.y()], hA=[p["hA_sn"]])

    T_hc4.set_dTdt_advective(source=T_hc3.y())
    T_hc4.dTdt_convective = T_hc3.dTdt_convective

    # core
    n.set_dndt(
        r=rho.y() + rho_ext,
        beta_eff=p["beta_t"],
        Lambda=p["Lam"],
        lam=p["lam"],
        C=[C1.y(), C2.y(), C3.y(), C4.y(), C5.y(), C6.y()]
    )
    C1.set_dcdt(n.y(), beta = p["beta"][0], Lambda = p["Lam"], lam = p["lam"][0], t_c = p["tau_c"], t_l = p["tau_l"], flow = True)
    C2.set_dcdt(n.y(), beta = p["beta"][1], Lambda = p["Lam"], lam = p["lam"][1], t_c = p["tau_c"], t_l = p["tau_l"], flow = True)
    C3.set_dcdt(n.y(), beta = p["beta"][2], Lambda = p["Lam"], lam = p["lam"][2], t_c = p["tau_c"], t_l = p["tau_l"], flow = True)
    C4.set_dcdt(n.y(), beta = p["beta"][3], Lambda = p["Lam"], lam = p["lam"][3], t_c = p["tau_c"], t_l = p["tau_l"], flow = True)
    C5.set_dcdt(n.y(), beta = p["beta"][4], Lambda = p["Lam"], lam = p["lam"][4], t_c = p["tau_c"], t_l = p["tau_l"], flow = True)
    C6.set_dcdt(n.y(), beta = p["beta"][5], Lambda = p["Lam"], lam = p["lam"][5], t_c = p["tau_c"], t_l = p["tau_l"], flow = True)

    T_cg.set_dTdt_convective(source=[T_cf1.y()], hA=[p["hA_fg"]])
    T_cg.set_dTdt_internal(source = [n.y()], k = [p["k_g"] * p["P"]])

    T_cf1.set_dTdt_advective(source=T_hf4.y(t - p["tau_hx_c"]))
    T_cf1.set_dTdt_convective(source=[T_cg.y()], hA=[p["k_1"] * p["hA_fg"]])
    T_cf1.set_dTdt_internal(source= [n.y()], k = [p["k_f1"] * p["P"]])

    T_cf2.set_dTdt_advective(source=T_cf1.y())
    T_cf2.dTdt_convective = T_cf1.dTdt_convective
    T_cf2.set_dTdt_internal(source=[n.y()], k=[p["k_f2"] * p["P"]])

    rho.set_drdt(
        sources=[T_cf1.dydt, T_cf2.dydt, T_cg.dydt],
        coeffs=[p["a_f"] / 2, p["a_f"] / 2, p["a_g"]]
    )

    return p, MSRE


# test msre model against simulink data from Singh et al. 2015 (https://doi.org/10.1016/j.anucene.2017.10.047)
def test_msre_model_1MW():

    # parameters
    P = 1.0
    inserted = 1.39e-4
    duration = 415.0
    t_ins = 2500.0

    # Load reference data
    df_simulink = pd.read_excel(f"{os.getcwd()}/tests/msre_data/simulink_msre_{int(P)}MW_U233_insertion.xlsx")
    i_trans = [i for i in range(len(df_simulink)) if df_simulink['time'][i] >= 2500]

    # generate solution using msrDynamics
    T = np.array(df_simulink['time'])
    p, model = _build_msre_model(power = P, times = T, t_ins = t_ins, pcm = inserted)
    _ = model.solve(T, max_delay = p['tau_l'], populate_nodes = True)
    n = model.nodes['n']

    # parse data for comparison
    i_insert = np.array([i for i in range(len(T)) if (T[i] > t_ins) and (T[i] < t_ins + duration)])
    ref_P = P*n.y_out[i_insert[0]-100]
    ref_P_simulink = df_simulink['Mux(4)'][i_trans[0]-100]*P
    P_dat = np.array([(k*P)-ref_P for k in n.y_out])
    dP_msrD = P_dat[i_insert]
    dP_simulink = df_simulink['Mux(4)'][i_insert]*P-ref_P_simulink

    # compare
    np.testing.assert_allclose(dP_msrD, dP_simulink, rtol=1e-1, atol = 1e-2)

def test_msre_model_5MW():

    # parameters
    P = 5.0
    inserted = 1.96e-4
    duration = 415.0
    t_ins = 2500.0

    # Load reference data
    df_simulink = pd.read_excel(f"{os.getcwd()}/tests/msre_data/simulink_msre_{int(P)}MW_U233_insertion.xlsx")
    i_trans = [i for i in range(len(df_simulink)) if df_simulink['time'][i] >= 2500]

    # generate solution using msrDynamics
    T = np.array(df_simulink['time'])
    p, model = _build_msre_model(power = P, times = T, t_ins = t_ins, pcm = inserted)
    _ = model.solve(T, max_delay = p['tau_l'], populate_nodes = True)
    n = model.nodes['n']

    # parse data for comparison
    i_insert = np.array([i for i in range(len(T)) if (T[i] > t_ins) and (T[i] < t_ins + duration)])
    ref_P = P*n.y_out[i_insert[0]-100]
    ref_P_simulink = df_simulink['Mux(4)'][i_trans[0]-100]*P
    P_dat = np.array([(k*P)-ref_P for k in n.y_out])
    dP_msrD = P_dat[i_insert]
    dP_simulink = df_simulink['Mux(4)'][i_insert]*P-ref_P_simulink

    # compare
    np.testing.assert_allclose(dP_msrD, dP_simulink, rtol=1e-1, atol = 5*1e-2)

def test_msre_model_8MW():

    # parameters
    P = 8.0
    inserted = 2.48e-4
    duration = 415.0
    t_ins = 2500.0

    # Load reference data
    df_simulink = pd.read_excel(f"{os.getcwd()}/tests/msre_data/simulink_msre_{int(P)}MW_U233_insertion.xlsx")
    i_trans = [i for i in range(len(df_simulink)) if df_simulink['time'][i] >= 2500]

    # generate solution using msrDynamics
    T = np.array(df_simulink['time'])
    p, model = _build_msre_model(power = P, times = T, t_ins = t_ins, pcm = inserted)
    _ = model.solve(T, max_delay = p['tau_l'], populate_nodes = True)
    n = model.nodes['n']

    # parse data for comparison
    i_insert = np.array([i for i in range(len(T)) if (T[i] > t_ins) and (T[i] < t_ins + duration)])
    ref_P = P*n.y_out[i_insert[0]-100]
    ref_P_simulink = df_simulink['Mux(4)'][i_trans[0]-100]*P
    P_dat = np.array([(k*P)-ref_P for k in n.y_out])
    dP_msrD = P_dat[i_insert]
    dP_simulink = df_simulink['Mux(4)'][i_insert]*P-ref_P_simulink

    # compare
    np.testing.assert_allclose(dP_msrD, dP_simulink, rtol=1e-1, atol = 8*1e-2)