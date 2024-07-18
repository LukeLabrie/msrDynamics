.. Project Name documentation master file, created by
   sphinx-quickstart on Fri Jun 18 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

msrDynamics
=========================================

*Documentation Under Construction*

``msrDynamics`` is an object-oriented API for constructing and solving reduced order thermal hydraulic systems, written 
with emulation of Simulink-style solvers for molten salt reactor (MSR) systems in mind 
(e.g. `Singh et al <https://www.sciencedirect.com/science/article/pii/S030645491730381X>`_), but can be extended to 
other fission and/or thermal hydraulic systems. The goal of this package is to streamline the implemetation of such 
nodal models for more complex systems, where direct handling of the equations can become cumbersome. 


Theory & Background
===================

The API is designed for building nodal systems where variables representing system components or masses are aggregated 
into nodes with associated properties, and which only interact with other nodes through these properties. The nodes 
therefore, are essentially representations of the state variables of a first-order system of differential equations, 
suitable for a numerical solver. 

A first-order system of ODEs (ordinary differential equations) takes the form

.. math::
   \dot{y} = f(t, y, \ldots)

For MSRs and other thermal hydraulic systems, the dynamics of the state at a given time, :math:`\dot{y(t)}`, depend on 
the state at time :math:`t`, :math:`y(t)`, but may also depend on the state at some previous time :math:`y(t-\tau)` 
where :math:`\tau = (\tau_1, \tau_2, ...)` are the respective delay times of the relevant state variables. The system to 
be solved is therefore a delay differential equation of the form 

.. math::

   \dot{y} = f(t, y(t), y(t-\tau_1), y(t-\tau_2), \ldots) 

``msrDynamics`` essentially provides tools to construct symbolic expressions representing a system of this form, which 
are then passed to `JiTCDDE <https://github.com/neurophysik/jitcdde>`_ as the backend solver. ``JiTCDDE`` integrates 
adaptively using a Bogacki-Shampine pair. The state and derivative of the system or stored at each integration step, 
which allows for past states to be evaluated using piece-wise cubic hermite inteprolation. This is the DDE integration 
method proposed by Shampine and Thompson [ST01]_.

Modeling Approach
=================

The model, as intended, essentially describes heat flow through the component masses, aggregating all spatial 
considerations into one-dimensional relationships governed by constant coefficients. Listed below is some relevant 
literature on this approach and its applications:

* `Ball, 1963 <https://digital.library.unt.edu/ark:/67531/metadc1201699/>`_
* `Ball and Kerlin, 1965 <https://www.osti.gov/biblio/4591881>`_
* `Kerlin, Ball and Steffy, 1971 <http://moltensalt.org/references/static/downloads/pdf/ORNL-TM-2571.pdf>`_
* `Singh et al., 2017 <https://doi.org/10.1016/j.anucene.2017.10.047>`_
* `Singh et al., 2018a <https://doi.org/10.1080/00295450.2017.1416879>`_
* `Singh et al., 2018b <https://doi.org/10.1016/j.anucene.2017.10.047>`_
* `Singh et al., 2020 <https://doi.org/10.1016/j.nucengdes.2019.110457>`_



Usage & API Description
=======================

The ``Node()`` object includes helper functions to define symbolic expressions representing convective and advective heat 
transfer, as well as generation from point-kinetics. User-defined dynamics are supported as well.  


Module Reference
================

.. automodule:: _msrDynamics

.. autoclass:: Node

.. autoclass:: System



=================
Backend
=================

`JiTCDDE <https://github.com/neurophysik/jitcdde>`_, a delay differential 
equation solver,


References
==========

.. [ST01] L.F. Shampine, S. Thompson: Solving DDEs in Matlab, Applied Numerical Mathematics 37, pp. 441–458 (2001), `10.1016/S0168-9274(00)00055-6 <http://dx.doi.org/10.1016/S0168-9274(00)00055-6>`_.

.. _JiTCDDE documentation: http://jitcdde.readthedocs.io