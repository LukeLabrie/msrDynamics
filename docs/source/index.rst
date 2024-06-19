.. Project Name documentation master file, created by
   sphinx-quickstart on Fri Jun 18 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

msrDynamics
=========================================

**Under Construction**

`msrDynamics` is an object-oriented API for constructing and solving reduced order thermal hydraulic systems, written 
with emulation of Simulink-style solvers for molten salt reactor (MSR) systems in mind 
(e.g. `Link Singh et al <https://www.sciencedirect.com/science/article/pii/S030645491730381X>`_), but can be extended to 
other fission and/or thermal hydraulic systems. The goal of this package is to streamline the implemetation of such 
nodal models for more complex systems, where direct handling of the equations can become cumbersome. 


Theoretical Background
======================

The API is designed for building nodal systems where variables representing system components or masses are aggregated 
into nodes with associated properties, and which only interact with other nodes through these properties. The nodes 
therefore, are essentially representations of the state variables of a first-order system of differential equations, 
suitable for a numerical solver. 

A first-order system of ODEs (ordinary differential equations) takes the form

.. math::
   \dot{y} = f(t, y)


Listed below is some relevant literature on this approach and its applications:
* `Ball, 1963 <https://digital.library.unt.edu/ark:/67531/metadc1201699/>`_
* `Ball and Kerlin, 1965 <https://www.osti.gov/biblio/4591881>`_
* `Kerlin, Ball and Steffy, 1971 <http://moltensalt.org/references/static/downloads/pdf/ORNL-TM-2571.pdf>`_
* `Singh et al., 2017 <https://doi.org/10.1016/j.anucene.2017.10.047>`_
* `Singh et al., 2018a <https://doi.org/10.1080/00295450.2017.1416879>`_
* `Singh et al., 2018b <https://doi.org/10.1016/j.anucene.2017.10.047>`_
* `Singh et al., 2020 <https://doi.org/10.1016/j.nucengdes.2019.110457>`_



Model & Usage
=============




=================
Backend
=================

to `Link JiTCDDE <https://github.com/neurophysik/jitcdde>`_, a delay
differential equation solver,

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   installation
   usage
   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
