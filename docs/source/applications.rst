Applications
============

Overview
--------

To showcase some of OpenSBLI's functionality, three applications are presented here.

1D wave propagation
-------------------

The one-dimensional wave equation is written as

.. math:: \frac{\partial^2 u}{\partial t^2} = c\frac{\partial^2 u}{\partial x^2},

where :math:`u`  is the dispacement and :math:`c`  is a known constant.

2D advection-diffusion of a Gaussian bump
-----------------------------------------



3D Taylor-Green vortex
----------------------

This problem solves the compressible Navier-Stokes equations in 3D. The equations governing conservation of mass, momentum and energy are respectively written as

.. math:: \frac{\partial \rho}{\partial t} + \frac{\partial}{\partial x_j}\left[\rho u_j\right] = 0,

.. math:: \frac{\partial \rho u_i}{\partial t} + \frac{\partial}{\partial x_j}\left[\rho u_i u_j + p\delta_{ij} - \tau_{ij}\right] = 0,

and

.. math:: \frac{\partial \rho E}{\partial t} + \frac{\partial}{\partial x_j}\left[\rho E u_j + u_j p + q_i - u_i\tau_{ij}\right] = 0.

The quantity :math:`\rho` is the fluid density, :math:`u` is the velocity, :math:`p` is the pressure field, :math:`E` is the total energy. The components of the stress tensor :math:`\tau` are given by

.. math:: \tau_{ij} = \frac{\mu}{\mathrm{Re}}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right) - \frac{2}{3}\frac{\mu}{\mathrm{Re}}\frac{\partial u_k}{\partial x_k}\delta_{ij},

where :math:`\delta_{ij}` is the Kronecker Delta function and :math:`\mathrm{Re}` is the Reynolds number. The components of the heat flux term :math:`q` are given by

.. math:: q_i = -\frac{\mu}{(\gamma-1)\ \mathrm{M}_\inf^2\ \mathrm{Pr}\ \mathrm{Re}}\frac{\partial T}{\partial x_i},

where :math:`T` is the temperature field, :math:`\gamma` is the ratio of specific heats, :math:`\mathrm{M}_\inf` is the free-stream Mach number, and :math:`\mathrm{Pr}` is the Prandtl number.
