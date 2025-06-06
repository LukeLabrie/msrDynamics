from ._msrDynamics import Node
from symengine import Function
import numpy as np

class PID_loop:
    """
    Implements a PID (Proportional-Integral-Derivative) control loop with options
    for customization, state tracking, and boundary constraints.

    Attributes:
        base_value (float): The base value added to the PID output.
        setpoint_node (Node): The node representing the setpoint.
        setpoint_value (float): The desired setpoint value.
        k_p (float): Proportional gain.
        k_i (float): Integral gain.
        k_d (float): Derivative gain.
        name (str): Identifier for the PID loop.
        n_args (int): Number of arguments required for the symbolic function.
        initial_value (float): Initial value for the PID controller.
        bound (tuple): Tuple specifying the lower and upper bounds for the output.
        clegg_integrator (bool): Enables resetting the integrator upon error sign change.
        min_reading (float): Minimum reading threshold to consider valid state input.
        output_sym (Function): Symbolic representation of the PID output.
        state (list): Stores state values over time.
        err (list): Stores error values over time.
        dt (list): Stores time step values over time.
        times (list): Stores time values over time.
        p_output (list): Stores proportional output values over time.
        i_output (list): Stores integral output values over time.
        d_output (list): Stores derivative output values over time.
        output (list): Stores PID output values over time.
        dedt (list): Stores derivative of error over time.
        integral (list): Stores cumulative integral of the error over time.
    """

    def __init__(self,
                 base_value: float,
                 setpoint_node: Node,
                 setpoint_value: float,
                 k_p: float = 0.0,
                 k_i: float = 0.0,
                 k_d: float = 0.0,
                 name: str = None,
                 n_args: int = 2,
                 initial_value: float = None,
                 bound: tuple = None,
                 clegg_integrator: bool = False,
                 min_reading: float = None,
                 store_all_output: bool = False,
                 ) -> None:
        """
        Initializes the PID_loop object.

        Args:
            base_value (float): The base value added to the PID output.
            setpoint_node (Node): The node representing the setpoint.
            setpoint_value (float): The desired setpoint value.
            k_p (float, optional): Proportional gain. Defaults to 0.0.
            k_i (float, optional): Integral gain. Defaults to 0.0.
            k_d (float, optional): Derivative gain. Defaults to 0.0.
            name (str, optional): Identifier for the PID loop. Defaults to None.
            n_args (int, optional): Number of arguments for the symbolic function. Defaults to 2.
            initial_value (float, optional): Initial value for the PID controller. Defaults to None.
            bound (tuple, optional): Tuple specifying lower and upper bounds for the output. Defaults to None.
            clegg_integrator (bool, optional): Enables resetting the integrator upon error sign change. Defaults to False.
            min_reading (float, optional): Minimum reading threshold for valid state input. Defaults to None.
        """

        self.base_value = base_value
        self.setpoint_node = setpoint_node
        self.setpoint_value = setpoint_value 
        self.k_p = k_p
        self.k_i = k_i
        self.k_d = k_d
        if name is None:
            self.name = f"pid_loop_{setpoint_node.name}_{setpoint_value}"
        else:
            self.name = name
        self.n_args = n_args
        if initial_value is None:
            self.initial_value = setpoint_node.y0
        else:
            self.initial_value = initial_value
        self.output_sym = Function(self.name)
        self._output_func = None
        self.cumsum = 0.0
        if store_all_output:
            self.times = []
            self.output = []
            self.state = []
            self.p_output = []
            self.i_output = []
            self.d_output = []
            self.err = []
            self.dt = []
            self.dedt = []
            self.integral = []
        else:
            self.times = None
            self.output = None
            self.state = None
            self.p_output = None
            self.i_output = None
            self.d_output = None
            self.err = None
            self.dt = None
            self.dedt = None
            self.integral = None

        self.err_prev = None
        self.bound = bound
        self.clegg_integrator = clegg_integrator
        self.de_prev = None
        self.min_reading = min_reading
        self.store_all_output = store_all_output
        self.t_prev = None
        self.out_prev = None

        def output_func(self):
            """
            Generates or returns the PID controller function.

            Returns:
                callable: A function implementing the PID control logic.
            """
            if self._output_func is None:
                def pid_func(y, state, t):
                    """
                    PID control logic for calculating the output.

                    Args:
                        y (float): Current output value.
                        state (float): Current state value.
                        t (float): Current time.

                    Returns:
                        float: PID control output value.
                    """

                    # unpack attributes for performance
                    t_prev = self.t_prev
                    min_reading = self.min_reading
                    setpoint_value = self.setpoint_value
                    err_prev = self.err_prev
                    ci = self.clegg_integrator
                    out_prev = self.out_prev
                    k_p = self.k_p
                    k_i = self.k_i
                    k_d = self.k_d
                    base_value = self.base_value
                    bound = self.bound
                    dt = t - t_prev if t_prev else t
                    if (min_reading) and (state < min_reading):
                        p_out, i_out, d_out, out, err, dedt = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                    else:
                        # p
                        err = state - setpoint_value
                        if err_prev and ci:
                            if np.sign(err_prev) != np.sign(err):
                                self.cumsum = 0.0
                        # i
                        self.cumsum += err*dt
                        # d
                        de = err - err_prev if err_prev is not None else 0.0
                        if dt == 0.0:
                            out = out_prev if out_prev else 0.0
                            return out
                        dedt = de / dt

                        p_out = k_p*err
                        i_out = k_i*self.cumsum
                        d_out = k_d*dedt
                        calc = p_out + d_out + i_out + base_value

                        if bound:
                            out = max(bound[0], min(calc, bound[1]))
                        else:
                            out = calc

                    # store inputs/outputs
                    self.err_prev = err
                    self.t_prev = t
                    self.out_prev = out
                    if self.store_all_output:
                        self.times.append(t)
                        self.output.append(out)
                        self.state.append(state)
                        self.p_output.append(p_out)
                        self.i_output.append(i_out)
                        self.d_output.append(d_out)
                        self.err.append(err)
                        self.dt.append(dt)
                        self.dedt.append(dedt)
                        self.integral.append(self.cumsum)
                    return out

                return pid_func
            
        self.output_func = output_func(self)
