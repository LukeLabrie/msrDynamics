import numpy as np
import pandas as pd

parameters = {
    # Domain
    "t0": 0.0,  # Start time (s)
    "tf": 1000.00,  # End time (s)
    "T": np.arange(0.0, 1000.00, 0.01),  # Time array (s)

    # Perturbations - SOURCE INSERTION
    "sourcedata": np.array([0, 0, 0]),  # Source data (neutron flux/s)
    "sourcetime": np.array([0, 50, 100]),  # Source time points (s)
    "source": pd.Series(np.array([0, 0, 0]), index=np.array([0, 50, 100])),

    # Perturbations - REACTIVITY INSERTION
    "simtime": 10,  # Simulation time (s)
    "reactdata": np.array([0, 5E-4]),  # Reactivity data (delta rho)
    "reacttime": np.array([0, 2500]),  # Reactivity time points (s)
    "react": pd.Series(np.array([0, 5E-4]), index=np.array([0, 2500])),
    "ts_max": 1e-1,  # Maximum timestep (s)

    # NEUTRONICS DATA
    "tau_l": 16.73,  # ORNL-TM-0728 (s)
    "tau_c": 8.46,  # ORNL-TM-0728 (s)
    "P": 8.0,  # Thermal Power in MW
    "n_frac0": 1.0,  # Initial fractional neutron density
    "Lam": 2.400E-04,  # Mean generation time (s)
    "lam": np.array([1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00]),  # Delayed neutron group decay constants
    "beta": np.array([0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023]),  # Delayed neutron fractions
    "beta_t": np.sum(np.array([0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023])),  # Total delayed neutron fraction
    "C0": None,  # Placeholder for dependent calculation
}

# Calculate dependent parameters
parameters["rho_0"] = parameters["beta_t"] - sum(
    np.divide(parameters["beta"], 1 + np.divide(1 - np.exp(-parameters["lam"] * parameters["tau_l"]), parameters["lam"] * parameters["tau_c"]))
)

parameters["C0"] = parameters["beta"] / parameters["Lam"] * (
    1.0 / (parameters["lam"] - (np.exp(-parameters["lam"] * parameters["tau_l"]) - 1.0) / parameters["tau_c"])
)

parameters.update({
    # Feedback coefficients
    "a_f": -8.71E-05,  # Fuel salt temperature-reactivity feedback coefficient
    "a_g": -6.66E-05,  # Graphite temperature-reactivity feedback coefficient

    # CORE HEAT TRANSFER PARAMETERS
    # Fuel Parameters
    "vdot_f": 7.5708E-02,  # Volumetric flow rate (m^3/s)
    "rho_f": 2.14647E+03,  # Density of fuel salt (kg/m^3)
    "W_f": 1.623879934566580e+02,  # Fuel flow rate (kg/s)
    "m_f": parameters["W_f"] * parameters["tau_c"],  # Fuel mass in core (kg)
    "nn_f": 2,  # Number of fuel nodes in core model
    "mn_f": parameters["m_f"] / 2,  # Fuel mass per node (kg)
    "scp_f": 1.9665E-3,  # Specific heat capacity of fuel salt (MJ/kg-C)

    # Core Upflow
    "v_g": 1.95386,  # Graphite volume (m^3)
    "rho_g": 1.860E3,  # Graphite density (kg/m^3)
    "m_g": 1.95386 * 1.860E3,  # Graphite mass (kg)
    "scp_g": 1.773E-3,  # Specific heat capacity of graphite (MW-s/kg-C)
    "hA_fg": 0.02 * 9 / 5,  # Heat transfer coefficient (MW/°C)
    "k_g": 0.07,  # Fraction of total power generated in graphite
    "k_f": 0.93,  # Fraction of heat generated in fuel
    "k_f1": 0.93 / 2,  # Fraction of power in fuel lump 1

    # Initial conditions
    "Tf_in": 6.3222E+02,  # Inlet temperature (°C)
    "T0_f2": 6.5727E+02,  # Fuel temperature lump 2 (°C)
    "T0_f1": 6.3222E+02 + (6.5727E+02 - 6.3222E+02) / 2,  # Fuel temperature lump 1 (°C)

    # Heat Exchanger
    "d_he": 16,  # Heat exchanger diameter (in)
    "h_he": 72,  # Heat exchanger height (in)
    "od_tube": 0.5,  # Coolant tube OD (in)
    "n_tube": 159,  # Number of coolant tubes

    # PRIMARY FLOW PARAMETERS
    "W_p": 1.623879934566580e+02,  # Fuel flow rate (kg/s)

    # SECONDARY FLOW PARAMETERS
    "vdot_s": 5.36265E-02,  # Coolant volumetric flow rate (m^3/s)
    "rho_s": 1.922e3,  # Coolant salt density (kg/m^3)

    # Radiator Parameters
    "od_rad": 0.01905,  # Outer diameter of tubes in radiator (m)
    "tube_wall_thick": 0.0018288,  # Tube wall thickness (m)

    # Initial conditions
    "T0_rp": 6.5727E+02,  # Primary side initial temperature (°C)

    # Pure time delays
    "tau_hx_c": 8.67,  # Delay from HX to core (s)
})

