import numpy as np

def get_params(power):

    params = {
        # Domain
        "t0": 0.0,  # Start time
        "tf": 5000.00,  # End time
        "T": np.arange(0.0, 5000.00, 0.01),  # Time mesh

        # Reactivity Insertion
        "inserted": 1.39e-4,  # Reactivity insertion for 1MW
        # "inserted": 1.96e-4,  # Reactivity insertion for 5MW
        # "inserted": 2.48e-4,  # Reactivity insertion for 8MW

        # Neutronics Data
        "tau_l": 16.73,  # Lumped parameter delay time (s)
        "tau_c": 8.46,  # Circulation time (s)
        "P": power,  # Reactor power (MW, example values are 0.1, 5, or 8)
        "n_frac0": 1,  # Initial fractional neutron density n/n0
        "Lam": 4.0E-04,  # Neutron generation time (s)
        "lam": np.array([1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00]),  # Decay constants (s^-1)
        "beta": np.array([0.00023, 0.00079, 0.00067, 0.00073, 0.00013, 0.00009]),  # Delayed neutron fractions
        "beta_t": np.sum(np.array([0.00023, 0.00079, 0.00067, 0.00073, 0.00013, 0.00009])),  # Total delayed neutron fraction
    }

    params.update({
        "rho_0": params["beta_t"] - np.sum(  # Reactivity change (stationary to circulating fuel)
            np.divide(
                params["beta"],
                1 + np.divide(
                    1 - np.exp(-params["lam"] * params["tau_l"]),
                    params["lam"] * params["tau_c"]
                )
            )
        ),
        "C0": params["beta"] / params["Lam"] * (  # Initial concentration of precursors
            1.0 / (
                params["lam"] - (np.exp(-params["lam"] * params["tau_l"]) - 1.0) / params["tau_c"]
            )
        )
    })

    params.update({
        # Feedback Coefficients
        "a_f": -11.034E-5,  # Fuel temperature feedback coefficient (pcm/K)
        "a_g": -05.814E-5,  # Graphite temperature feedback coefficient (pcm/K)

        # Core Heat Transfer Parameters
        "vdot_f": 7.5708E-02,  # Core fuel flow rate (m^3/s)
        "rho_f": 2.14647E+03,  # Fuel density (kg/m^3)
        "W_f": 1.623879934566580e+02,  # Core fuel mass flow rate (kg/s)
    })

    params.update({
        "m_f": params["W_f"] * params["tau_c"],  # Total fuel mass in core (kg)
        "nn_f": 2,  # Number of fuel channels
    })

    params.update({
        "mn_f": params["m_f"] / params["nn_f"],  # Mass per fuel channel (kg)
        "scp_f": 1.9665E-3,  # Specific heat capacity of fuel (kJ/kg-K)

        # Core Upflow
        "v_g": 1.95386,  # Graphite core flow velocity (m/s)
        "rho_g": 1.860E3,  # Graphite density (kg/m^3)
    })

    params.update({
        "m_g": params["v_g"] * params["rho_g"],  # Graphite core mass (kg)
        "scp_g": 1.773E-3,  # Specific heat capacity of graphite (kJ/kg-K)
    })

    params.update({
        "mcp_g1": params["m_g"] * params["scp_g"],  # Heat capacity of graphite (kJ/K)
        "mcp_f1": params["mn_f"] * params["scp_f"],  # Heat capacity of fuel in channel 1 (kJ/K)
        "mcp_f2": params["mn_f"] * params["scp_f"],  # Heat capacity of fuel in channel 2 (kJ/K)
        "hA_fg": 0.02 * 9 / 5,  # Heat transfer coefficient fuel-to-graphite (kJ/s-K)
        "k_g": 0.07,  # Thermal conductivity of graphite (kW/m-K)
        "k_1": 0.5,  # Thermal conductivity parameter 1
        "k_2": 0.5,  # Thermal conductivity parameter 2
        "k_f": 0.93,  # Thermal conductivity of fuel (kW/m-K)
    })

    params.update({
        "k_f1": params["k_f"] / params["nn_f"],  # Thermal conductivity per fuel channel 1
        "k_f2": params["k_f"] / params["nn_f"],  # Thermal conductivity per fuel channel 2

        # Heat Exchanger Parameters
        "d_he": 16,  # Diameter of heat exchanger (m)
        "h_he": 72,  # Height of heat exchanger (m)
        "od_tube": 0.5,  # Outer diameter of tubes (m)
    })

    params.update({
        "id_tube": params["od_tube"] - 2 * 0.042,  # Inner diameter of tubes (m)
        "n_tube": 159,  # Number of tubes
        "a_tube": 254 * 144,  # Total tube surface area (m^2)
    })

    params.update({
        "l_tube": params["a_tube"] / params["n_tube"] / (np.pi * params["od_tube"]),  # Length of tubes (m)
    })

    params.update({
        "v_tube": params["n_tube"] * np.pi * (params["od_tube"] / 2) ** 2 * params["l_tube"],  # Volume of tubes (m^3)
    })

    params.update({
        "v_cool": params["n_tube"] * np.pi * (params["id_tube"] / 2) ** 2 * params["l_tube"],  # Coolant volume (m^3)
        "v_he": (params["d_he"] / 2) ** 2 * np.pi * params["h_he"],  # Total heat exchanger volume (m^3)
    })

    params.update({
        "v_he_fuel": params["v_he"] - params["v_tube"],  # Fuel volume in heat exchanger (m^3)
        "in_m": 1.63871e-5,  # Inlet mass fraction
        "W_p": params["W_f"],  # Pumped flow rate (kg/s)
    })

    params.update({
        "m_p": params["v_he_fuel"] * params["in_m"] * params["rho_f"],  # Heat exchanger fuel mass (kg)
        "nn_p": 4,  # Number of fuel loops
    })

    params.update({
        "mn_p": params["m_p"] / params["nn_p"],  # Mass per fuel loop (kg)
        "cp_p": params["scp_f"],  # Specific heat capacity of pumped fuel (kJ/kg-K)
    })


    params.update({
        "vdot_s": 5.36265E-02,
        "rho_s": 1.922e3,
        "W_s": 1.005793369810108e+02,
    })

    params.update({
        "m_s": params["v_cool"] * params["in_m"] * params["rho_s"],
        "nn_s": 4,

    })

    params.update({
        "mn_s": params["m_s"] / params["nn_s"],
        "scp_s": 2.39E-3,
        "A_phe": 2.359E+01,
        "ha_p": 6.480E-01,
        "ha_s": 3.060E-01,
        "mcp_pn": params["mn_p"] * params["cp_p"],
    })

    params.update({
        "hA_pn": params["ha_p"] / params["nn_s"],
        "nn_t": 2,
        "rho_tube": 8.7745E+03,
    })

    params.update({
        "m_tn": (params["v_tube"] - params["v_cool"]) * params["in_m"] * params["rho_tube"] / params["nn_t"],
        "scp_t": 5.778E-04,
    })

    params.update({
        "mcp_tn": params["m_tn"] * params["scp_t"],
        "mcp_sn": params["mn_s"] * params["scp_s"],
        "hA_sn": params["ha_s"] / params["nn_s"]
    })

    # Initial conditions
    params.update({
        "Tf_in": 6.3222E+02,
        "T0_f2": 6.5727E+02,
    })
    params.update({
        "T0_f1": params["Tf_in"] + (params["T0_f2"] - params["Tf_in"]) / 2,
    })

    params.update({
        "T0_g1": params["T0_f1"] + ((params["k_g"] * params["P"] / params["hA_fg"])),
        "Tp_in": params["T0_f2"],
        "T0_p4": params["Tf_in"],
    })

    params.update({
        "T0_p1": params["Tp_in"] - (params["Tp_in"] - params["Tf_in"]) / 4,
        "T0_p2": params["Tp_in"] - 2 * (params["Tp_in"] - params["Tf_in"]) / 4,
        "T0_p3": params["Tp_in"] - 3 * (params["Tp_in"] - params["Tf_in"]) / 4,
        "Ts_in": 5.4611E+02,
        "T0_s4": 5.7939E+02,
    })

    params.update({
        "T0_s1": params["Ts_in"] + (params["T0_s4"] - params["Ts_in"]) / params["nn_s"],
        "T0_s2": params["Ts_in"] + 2 * (params["T0_s4"] - params["Ts_in"]) / params["nn_s"],
        "T0_s3": params["Ts_in"] + 3 * (params["T0_s4"] - params["Ts_in"]) / params["nn_s"],
    })

    params.update({
        "T0_t1": (params["T0_p1"] * params["hA_pn"] + params["T0_s3"] * params["hA_sn"]) / (params["hA_pn"] + params["hA_sn"]),
        "T0_t2": (params["T0_p3"] * params["hA_pn"] + params["T0_s1"] * params["hA_sn"]) / (params["hA_pn"] + params["hA_sn"]),
    })
    params.update({
        # Radiator Parameters
        "Trp_in": params["T0_s4"],  # Radiator primary inlet temperature (K)
        "T0_rp": params["Ts_in"],  # Radiator primary outlet temperature (K)
        "Trs_in": 37.78,  # Radiator secondary inlet temperature (K)
        "T0_rs": 148.9,  # Radiator secondary outlet temperature (K)
        "od_rad": 0.01905,  # Outer diameter of radiator tubes (m)
        "tube_wall_thick": 0.0018288,  # Tube wall thickness (m)
    })

    params.update({
        "id_rad": params["od_rad"] - 2 * params["tube_wall_thick"],  # Inner diameter of radiator tubes (m)
        "n_rtubes": 120,  # Number of radiator tubes
        "l_rtube": 9.144,  # Length of radiator tubes (m)
    })

    params.update({
        "v_rp": np.pi * (params["id_rad"] / 2) ** 2 * params["l_rtube"] * params["n_rtubes"],  # Radiator primary volume (m^3)
        "n_tpr": 12,  # Number of tubes per row
        "n_row": 10,  # Number of rows of tubes
        "tube_space": 0.0381,  # Spacing between tubes (m)
    })

    params.update({
        "v_rs": (params["n_row"] * params["od_rad"] + (params["n_row"] - 1) * params["tube_space"]) * \
                (params["n_tpr"] * params["od_rad"] + (params["n_tpr"] - 1) * params["tube_space"]) * params["l_rtube"],  # Radiator secondary volume (m^3)
        "W_rp": params["W_s"],  # Radiator primary flow rate (kg/s)
        "m_rp": params["v_rp"] * params["rho_s"],  # Radiator primary mass (kg)
        "nn_rp": 1,  # Number of primary loops
    })

    params.update({
        "mn_rp": params["m_rp"] / params["nn_rp"],  # Radiator primary mass per loop (kg)
        "cp_rp": params["scp_s"],  # Specific heat capacity of radiator primary fluid (kJ/kg-K)
        "vdot_rs": 94.389,  # Radiator secondary flow rate (m^3/s)
        "rho_rs": 1.1237,  # Density of radiator secondary fluid (kg/m^3)
    })

    params.update({
        "W_rs": params["vdot_rs"] * params["rho_rs"],  # Radiator secondary flow rate (kg/s)
        "m_rs": params["v_rs"] * params["rho_rs"],  # Radiator secondary mass (kg)
        "nn_rs": 1,  # Number of secondary loops
    })

    params.update({
        "mn_rs": params["m_rs"] / params["nn_rs"],  # Radiator secondary mass per loop (kg)
        "scp_rs": 1.0085E-3,  # Specific heat capacity of radiator secondary fluid (kJ/kg-K)
        "A_rad": 6.503E1,  # Radiator heat exchange area (m^2)
    })

    params.update({
        "mcp_rpn" : params["mn_rp"] * params["cp_rp"],  # Radiator primary heat capacity (kJ/K)
    })

    params.update({
        "h_roverall": params["P"] / params["A_rad"] / ((params["T0_rp"] + params["Trp_in"]) / 2 - (params["T0_rs"] + params["Trs_in"]) / 2),  # Overall heat transfer coefficient (kJ/s-K)
        "mcp_rpn": params["mn_rp"] * params["cp_rp"],  # Heat capacity of radiator primary fluid (kJ/K)
        "mcp_rsn": params["mn_rs"] * params["scp_rs"],  # Heat capacity of radiator secondary fluid (kJ/K)
    })

    params.update({
        "hA_rpn": params["h_roverall"] * params["A_rad"] / params["nn_rs"],  # Heat transfer radiator primary (kJ/s-K)
        "hA_rsn": params["h_roverall"] * params["A_rad"] / params["nn_rs"],  # Heat transfer radiator secondary (kJ/s-K)
    })
    params.update({
        # Pure Time Delays Between Components
        "tau_hx_c": 8.67,  # Heat exchanger core-to-secondary delay (s)
        "tau_c_hx": 3.77,  # Heat exchanger secondary-to-core delay (s)
        "tau_hx_r": 4.71,  # Heat exchanger to radiator delay (s)
        "tau_r_hx": 8.24  # Radiator to heat exchanger delay (s)
    })

    return params

