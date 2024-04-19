# msrDynamics

`msrDynamics` is an object-oriented API to [JiTCDDE](https://github.com/neurophysik/jitcdde), 
a delay differential equation solver, written with emulation of simulink-style solvers for
molten salt reactor (MSR) systems in mind (see [Singh et al](https://www.sciencedirect.com/science/article/pii/S030645491730381X)),
but can be extended to other fission and/or thermal hydraulic systems. The goal of this package is to streamline the implemetation of such nodal models for more complex systems, where direct handling of the equations can become cumbersome. 

Comparison of results generated with `msrDynamics` against similar work by [Singh et al.](https://doi.org/10.1016/j.anucene.2017.10.047)
 ![](figures/step_insertion_5MW.png) 


## Installation

Clone this repository. Then, within your preferred python environment, install the project with pip by providing the path to the project folder. 
```
pip install </path/to/msrDynamics>
```
The project is not yet available on PyPI. 

## Methodology 

The API creates nodal systems, of the form discussed in the examples below, whereby variables representing system properties or masses are aggregated into nodes with associated properties, and which only interact with other nodes through these properties. The system can then be described as a first-order system of differential equations, suitable for a numerical solver. The nodes therefore, are essentially representations of the state variables of the system. 

## Usage

Nodes are represented by the `Node()` object which stores properties associated with the node, and contains helper methods to define certian dynamics like convective and advective heat transfer as well as neutron kinetics. For example, fuel flow through a core in direct contact with a moderator material (e.g. graphite), could be set up as follows. 

```python
import paramerers
import numpy as np

# instantiate system 
msr = System()

# define nodes with mass m, scpecific heat capacity scp, and mass flow W
f1 = Node(m = m_f1, scp = scp_f, W = W_f)
f2 = Node(m = m_f2, scp = scp_f, W = W_f)
g  = Node(m = m_g, scp = scp_g)

# add nodes to system 
msr.add_nodes([f1,f2,g])

# define dynamics 
f1.set_dTdt_advective(source = f_in)
f1.set_dTdt_convective(source = g.y(), hA = [hA_fg])

f2.set_dTdt_advective(source = f_1.y())
f2.set_dTdt_convective(source = g.y(), hA = [hA_fg])

g.set_dTdt_convective(source = [f1.y(), f2.y()], hA = [hA_fg, hA_fg])

# solve 
T = np.arange(0,100,0.01)
msr.solve(T)
```

Note, for any system, a `System()` object is required for proper handling of the global indexing required for the [JiTCDDE](https://github.com/neurophysik/jitcdde) backend. **Nodes need to be added to the system object before dynamics are defined**. This is because certain global system information is required in order to index the variables properly for the backend. 

Nodes can also represent state variables associated with neutron kinetics, like neutron concentration $n(t)$ delayed neutorn precursor concentrations $C_i(t)$, and therefore does not need to be associated with a thermal mass. The helper methods take other nodes to which a given node is coupled, along with relevant system parameters, as arguments, and sets up the symbolic expressions representing the equations govenring the dynamics of the node. These symbolic expressions are the input to the [JiTCDDE](https://github.com/neurophysik/jitcdde) backend. Alternatively, if the user wishes to circumvent the helper methods, or would like to define other dynamics, a node's dynamics can be set directly with a user-defined symbolic expression through the `Node.dydt` attribute, e.g.

```python
# instantiate system 
f = System()

# define nodes
x = Node(m = m_x)
y = Node(m = m_y)

# add nodes to system
f.add_nodes([x,y])

# define dynamics 
x.dydt = x.y() - y.y()
y.dydt = - x.y() + y.y() 
```

## Examples

See the notebooks below for more detailed examples of usage, as well as comparison to experimental data.

- [Simple Reactor](./examples/toy_reactor)
- [Aircraft Reactor Experiment](./examples/are)
- [Molten Salt Reactor Experiment](./examples/msre)
