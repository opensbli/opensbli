#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

from .equations import EinsteinTerm
from .scheme import *
from .kernel import *


class TemporalDiscretisation(object):

    """ Perform a temporal discretisation of the equations on the numerical grid of solution points. """

    def __init__(self, temporal_scheme, grid, const_dt, spatial_discretisation):
        """ Formulate the time discretisation scheme as a series of computational kernels.

        :arg temporal_scheme: The time discretisation scheme.
        :arg grid: The numerical Grid of solution points.
        :arg bool const_dt: True if the time-step is constant, and False otherwise.
        :arg spatial_discretisation: The object that performs the spatial discretisation.
        :returns: None
        """

        # Time-stepping scheme
        self.scheme = temporal_scheme

        # Constant or variable time-step
        if const_dt:
            dt = EinsteinTerm('deltat')
            dt.is_constant = True
            dt.is_commutative = True
        else:
            raise NotImplementedError("Varying delta t is not implemented in the code.")

        # Coefficients required for the time-stepping scheme, from the relevant Butcher tableau.
        self.nstages = temporal_scheme.order
        if isinstance(temporal_scheme, ForwardEuler) and self.nstages == 1:
            self.coeff = None
        elif isinstance(temporal_scheme, RungeKutta) and self.nstages == 3:
            self.coeff = self.scheme.get_coefficients()
        else:
            raise ValueError("Only first-order Forward or third-order Runge-Kutta temporal discretisation schemes are allowed.")

        # Start computations: Any computations at the start of the time-step. Generally these are the 'save' equations.
        self.start_computations = []

        # Computations: The main computations to perform during each time-step.
        self.computations = []

        # The prognostic variables (e.g. rho, rhou0, rhou1, rhoE in the Navier-Stokes equations) that should be updated.
        self.prognostic_variables = []

        # End computations: As of now, no end computations, any end computations like, diagnostics should be updated
        # to this attribute for code generation
        self.end_computations = None

        # The residual arrays that contain the change in the RHS of each equation.
        out = []
        for residual in spatial_discretisation.residual_arrays:
            out.append(self.time_derivative(residual.keys()[0].args[0], dt, residual[residual.keys()[0]], grid))
            self.prognostic_variables.append(residual.keys()[0].args[0])

        # Formulate each step of the time-stepping scheme here as a computational Kernel.
        range_of_evaluation = [tuple([0, s]) for i, s in enumerate(grid.shape)] # Grid point index 0 to nx (or ny or nz)
        if self.nstages != 1:
            # The 'save' equations.
            start = [o[-1] for o in out]
            range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
            self.start_computations.append(Kernel(start, range_of_evaluation, "Save equations", grid))

            # The 'update' equations of the variables at time 't + k', where k is the Runge-Kutta loop iteration.
            equations = [o[0] for o in out]
            self.computations.append(Kernel(equations, range_of_evaluation, "RK new (subloop) update", grid))
            equations = [o[1] for o in out]
            self.computations.append(Kernel(equations, range_of_evaluation, "RK old update", grid))
        else:
            self.start_computations = None
            self.computations.append(Kernel(out, range_of_evaluation, "Euler update", grid))

        # Copy the LHS vectors and scalars in spatial discretisation to prognostic_classified
        self.prognostic_classified = spatial_discretisation.lhs_vectors

        return

    def time_derivative(self, function, dt, residual, grid):
        """ Return the equation(s) used to advance the model equations forward in time.

        :arg function: The function that the LHS time derivative operates on.
        :arg dt: The time-step.
        :arg residual: The residual (i.e. the change in the RHS of the equation).
        :arg grid: The grid of solution points.
        :returns: An equation, or list of equations (if there is more than one step involved in the temporal scheme), used to advance the model equations forward in time.
        :rtype: sympy.Eq, or list of sympy.Eq
        """

        if self.nstages == 1:
            equation = Eq(function, function + dt*residual, evaluate=False)
        elif self.nstages == 3:
            old = grid.work_array('%s_old' % function.base)
            equation_function = Eq(function, old + dt*self.scheme.new*residual, evaluate=False)
            equation_old = Eq(old, old + dt*self.scheme.old*residual, evaluate=False)
            save_equation = Eq(old, function)
            equation = [equation_function, equation_old, save_equation]
        return equation


class RungeKutta(Scheme):

    """ Runge-Kutta time-stepping scheme. """

    def __init__(self, order):
        """ Set up the Runge-Kutta stages and the coefficients.

        :arg int order: The order of accuracy of the scheme.
        """

        Scheme.__init__(self, "RungeKutta", order)

        self.stage = Symbol('stage', integer=True)
        self.old = IndexedBase('rkold')
        self.old.is_grid = False
        self.old.is_constant = True
        self.old.ranges = order
        self.new = IndexedBase('rknew')
        self.new.is_grid = False
        self.new.is_constant = True
        self.new.ranges = order
        self.old = self.old[self.stage]
        self.new = self.new[self.stage]

        return

    def get_coefficients(self):
        """ Return the coefficients of the Runge-Kutta update equations.

        :returns: A dictionary of (update_equation, coefficients) pairs.
        :rtype: dict
        """

        if self.order == 3:
            coeffs = {}
            coeffs[self.old.base] = [Rational(1.0,4.0), Rational(3.0,20), Rational(3.0,5.0)]
            coeffs[self.new.base] = [Rational(2,3), Rational(5,12), Rational(3,5)]
        return coeffs


class ForwardEuler(Scheme):

    """ First-order forward/explicit Euler time-stepping scheme. """

    def __init__(self):
        """ Set up the forward Euler scheme. """

        Scheme.__init__(self, "ForwardEuler", 1)

        return
