from ._msrDynamics import Node
from symengine import Function
import numpy as np

class PID_loop:

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
                 ) -> None:

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
        self.p_output = []
        self.i_output = []
        self.d_output = []
        self.output = []
        self.err = []
        self.dt = []
        self.times = []
        self.err_prev = None
        self.dedt = []
        self.state = []
        self.bound = bound
        self.clegg_integrator = clegg_integrator
        self.de_prev = None
        self.min_reading = min_reading
        self.integral = []
        
    @property
    def output_func(self):
        if self._output_func is None:
            def pid_func(y, state, t):

                dt = t - self.times[-1] if self.times else t
                if (self.min_reading) and (state < self.min_reading):
                    p_out, i_out, d_out, out, err, dedt = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                else:
                    # p
                    err = state - self.setpoint_value
                    if self.err_prev and self.clegg_integrator:
                        if np.sign(self.err_prev) != np.sign(err):
                            self.cumsum = 0.0
                    # i
                    self.cumsum += err*dt
                    # d
                    de = err - self.err_prev if self.err_prev is not None else 0.0
                    if dt == 0.0:
                        out = self.output[-1] if self.output else 0.0
                        return out
                    dedt = de / dt

                    p_out = self.k_p*err
                    i_out = self.k_i*self.cumsum
                    d_out = self.k_d*dedt
                    calc = p_out + d_out + i_out + self.base_value

                    if self.bound:
                        out = max(self.bound[0], min(calc,self.bound[1]))
                    else:
                        out = calc

                # store inputs/outputs
                self.err_prev = err
                self.state.append(state)
                self.times.append(t)
                self.p_output.append(p_out)
                self.i_output.append(i_out)
                self.d_output.append(d_out)
                self.output.append(out)
                self.err.append(err)
                self.dt.append(dt)
                self.dedt.append(dedt)
                self.integral.append(self.cumsum)
                return out

            return pid_func
        else:
            return self._output_func
            
            
    @output_func.setter
    def output_func(self, custom_output_func):
        self._output_func = custom_output_func



    