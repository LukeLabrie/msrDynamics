from ._msre_parameters import parameters as p 
import numpy as np
from jitcdde import t
from msrDynamics import Node, System
import os
import pandas as pd

# msre model builder using msrDynamics
def _build_msre_model():

    # instantiate system object
    MSRE = System()

    # radiator
    T_out_rc = Node(m=p["mn_rp"], scp=p["mcp_rpn"] / p["mn_rp"], W=p["W_rp"], y0=p["T0_rp"])
    T_out_air = Node(m=p["mn_rs"], scp=p["mcp_rsn"] / p["mn_rs"], W=p["W_rs"], y0=p["T0_rs"])

    # heat exchanger
    T_hf1 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p1"])
    T_hf2 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p2"])
    T_hf3 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p3"])
    T_hf4 = Node(m=p["mn_p"], scp=p["mcp_pn"] / p["mn_p"], W=p["W_p"], y0=p["T0_p4"])
    T_ht1 = Node(m=p["m_tn"], scp=p["scp_t"], y0=p["T0_t1"])
    T_ht2 = Node(m=p["m_tn"], scp=p["scp_t"], y0=p["T0_t2"])
    T_hc1 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s1"])
    T_hc2 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s2"])
    T_hc3 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s3"])
    T_hc4 = Node(m=p["mn_s"], scp=p["mcp_sn"] / p["mn_s"], W=p["W_s"], y0=p["T0_s4"])

    # core
    n = Node(y0=p["n_frac0"])
    C1 = Node(y0=p["C0"][0])
    C2 = Node(y0=p["C0"][1])
    C3 = Node(y0=p["C0"][2])
    C4 = Node(y0=p["C0"][3])
    C5 = Node(y0=p["C0"][4])
    C6 = Node(y0=p["C0"][5])
    rho = Node(y0=0.0)

    # add reactivity input
    t_ins = 500
    inserted = 5e-4

    def rho_insert(t):
        if t < t_ins:
            return 0.0
        else:
            return inserted

    rho_ext = MSRE.add_input(rho_insert, p["T"])

    T_cg = Node(m=p["mcp_g1"] / p["scp_g"], scp=p["scp_g"], y0=p["T0_g1"])
    T_cf1 = Node(m=p["mn_f"], scp=p["scp_f"], W=p["W_f"], y0=p["T0_f1"])
    T_cf2 = Node(m=p["mn_f"], scp=p["scp_f"], W=p["W_f"], y0=p["T0_f2"])

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

    return MSRE


# test msre model against simulink data from Singh et al. 2015 (https://doi.org/10.1016/j.anucene.2017.10.047)
def test_msre_model():

    model = _build_msre_model()
    sol = model.solve(0.0, 1000.0, 0.01, max_delay = p['tau_l'])

    



    # Load reference data
    data_dir = 'msre_data'
    reference_data = pd.read_csv(os.path.join(data_dir, 'reference_data.csv'))

    # Extract model results
    model_results = MSRE.get_results()

    # Compare model results to reference data
    for column in reference_data.columns:
        np.testing.assert_allclose(model_results[column], reference_data[column], rtol=1e-5, atol=1e-8)

    print("Regression test passed.")
