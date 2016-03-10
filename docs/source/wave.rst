Application: 1D wave propagation
================================

Equations
---------

This test case solves the numerical solution of the one-dimensional wave equation, written as

.. math:: \frac{\partial^2 u}{\partial t^2} = c\frac{\partial^2 u}{\partial x^2},

where :math:`u`  is the dispacement (to be found) and :math:`c`  is a known constant (set to 0.25 in this simulation).

Simulation setup
----------------

A domain of length :math:`0 \leq x \leq 1` m is considered, with grid spacing :math:`dx` = 0.1 m, and periodic boundaries. A fourth-order accurate central differencing scheme is used to spatially discretise the domain, and a third-order Runge-Kutta timestepping scheme is used to march the equation forward in time.

The initial condition is defined by

.. math:: u(x, t=0) = \sin(2\pix),

The simulation was run with a timestep of :math:`dt` = 0.01 s until time :math:`t` = 1 s.

Results
-------

