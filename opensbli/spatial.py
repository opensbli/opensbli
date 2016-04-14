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

import os
import subprocess
from sympy import *
from sympy.parsing.sympy_parser import *

import logging
LOG = logging.getLogger(__name__)

from .equations import *
from .utils import *
from .kernel import *
from .scheme import *


class Central(Scheme):

    """ Spatial discretisation scheme using central differences. """

    def __init__(self, order):
        """ Set up the scheme.

        :arg int order: The order of accuracy of the scheme.
        """
        Scheme.__init__(self, "Central", order)
        self.points = list(i for i in range(-order/2, order/2+1))
        return


class SpatialDerivative(object):

    """ The spatial derivatives of an arbitrary function 'F'
    on the numerical grid with the provided spatial scheme.

    For a wall boundary condition this will have a dependency on the grid range. """

    def __init__(self, spatial_scheme, grid, max_order):
        """ Initialise the spatial derivative, which gives the equations
        of spatial Derivatives for the various combinations of the spatial scheme and order of accuracy.

        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :arg int max_order: The maximum order of the derivative in the function.
        :returns: None
        """
        
        self.stencil = self.create_stencil(spatial_scheme, grid)
        self.points = spatial_scheme.points

        self.derivatives = []
        self.derivative_direction = grid.indices
        self.index_mapping = grid.mapped_indices
        self.deltas = grid.deltas
        base = IndexedBase('f',shape = grid.shape)
        base.is_grid = True; base.is_constant = False
        fn = base[grid.indices]
        self.fn = fn

        self.create_derivative_formulas(fn, max_order, grid)
        return

    def create_stencil(self, spatial_scheme, grid):
        """ Create the computational stencil referenced by the grid indices.

        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :returns: A list of solution points indexed on the grid describing the computational stencil.
        :rtype: list
        """

        stencil = [[] for dim in grid.shape]

        for dim, val in enumerate(grid.shape):
            if isinstance(spatial_scheme, Central):
                points = list(i for i in range(-spatial_scheme.order/2, spatial_scheme.order/2+1)) # The local indices of each point in the stencil (in dimension 'dim').
                grid.halos.append((-spatial_scheme.order/2, spatial_scheme.order/2)) # The range of the indices of the stencil at the boundary which includes the halo points.
            else:
                raise NotImplementedError("Only central difference schemes are supported.")
            stencil[dim] = [grid.indices[dim] + i for i in points] # The indices in the grid offset by the local indices of each stencil point (e.g. i-1, i, i+1)
        return stencil

    def create_derivative_formulas(self, fn, max_order, grid):
        """ Create the formulas for the derivatives of the given function,
        based on the stencil pattern provided.

        :arg fn: The function whose differential needs computing.
        :arg int max_order: The maximum order of the derivative in the function.
        :arg grid: The numerical grid of solution points.
        :returns: None
        """

        derivatives = []
        derivatives += [fn] # FIXME: Later change this to interpolation
        derivative_formula = []
        derivative_formula += [fn]
        kernels = []
        kernels += [fn]

        # Create derivatives for each order, up to the maximum order of derivative that is present in the equation.
        for order in range(1, max_order+1):
            shape = tuple([len(grid.indices) for ind in range(order)])
            array = MutableDenseNDimArray.zeros(*shape)
            fdarray = MutableDenseNDimArray.zeros(*shape)
            derivative_kernel = MutableDenseNDimArray.zeros(*shape)

            for ind in np.ndindex(*array.shape):
                # Arguments to the derivative in terms of grid indices
                derivative_args = [grid.indices[i] for i in ind]
                # The derivative kernel's name
                name = [str(arg) for arg in ind]
                name = "[%d][%s]" % (order, ','.join(name))
                derivative_kernel[ind] = Symbol(name)

                # Find the finite difference formula
                array[ind] = fn.diff(*derivative_args)
                if order == 1 or len(set(derivative_args)) == 1:
                    fdarray[ind] = as_finite_diff(array[ind], self.stencil[ind[0]])*pow(grid.deltas[ind[0]],-order)
                else:
                    new_derivative = array[ind].subs(derivatives[order-1][ind[:-1]], derivative_formula[order-1][ind[:-1]])
                    fdarray[ind] = as_finite_diff(new_derivative, self.stencil[ind[-1]], wrt=grid.indices[ind[-1]])*pow(grid.deltas[ind[-1]],-1)

            derivatives.append(array)
            derivative_formula.append(fdarray)
            kernels.append(derivative_kernel)

        self.derivatives = derivatives
        self.derivative_formula = derivative_formula
        self.derivative_kernel = kernels
        return

    def get_derivative_formula(self, derivative):
        """
        This returns the formula for the derivative function
        """
        order = len(derivative.args[1:])
        indices = list(derivative.args[1:])
        if order == 1 or len(set(indices)) == 1:
            wrt = indices[0]
            gridpoints = [wrt + i for i in self.points]
            formula = apply_finite_diff(order, gridpoints, [derivative.expr.subs({wrt: x}) for x in gridpoints], wrt)
            d1 = self.derivative_direction.index(self.index_mapping[derivative.args[1]])
            delta = self.deltas[d1]
            formula = formula*pow(delta, -order)
        elif order == 2:
            # Do the derivative of derivative
            raise NotImplementedError("Derivatives of order == 2  Mixed is not implemented")
        else:
            raise NotImplementedError("Derivatives of order > 2 are not implemented")
        return formula
    def get_derivative(self, derivative):
        """ Return a tuple to which the derivative formula exists in
        the already-evaluated derivatives.

        :arg derivative: The derivative you want to get the formula for.
        """
        order = len(derivative.args[1:])
        indices = []
        for arg in derivative.args[1:]:
            indices = indices + [self.derivative_direction.index(self.index_mapping[arg])]
        general_formula = []
        subevals = []
        requires = []
        if order == 1 or len(set(indices)) == 1:
            general_formula += [order,tuple(indices)]
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += list(derivative.args[0].atoms(Indexed))
        else:
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += [derivative.args[0]]
            general_formula += [order-1, tuple([indices[-1]])]
            requires += [Derivative(derivative.args[0],*derivative.args[1:-1])]

        return general_formula, subevals, requires


class SpatialDiscretisation(object):

    """ The spatial discretisation using the provided scheme on the provided grid. """

    def __init__(self, expanded_equations, expanded_formulas, grid, spatial_scheme):
        """ Perform the spatial discretisation.

        :arg list expanded_equations: A list of the equations expanded with respect to the Einstein indices.
        :arg list expanded_formulas: A list of the formulas expanded with respect to the Einstein indices.
        :arg grid: The numerical grid of solution points.
        :arg spatial_scheme: The spatial scheme used to perform the spatial discretisation.
        :returns: None
        """

        all_equations = flatten(expanded_equations)
        all_formulas = flatten(expanded_formulas)
        all_formulas = get_used_formulas(all_formulas, all_equations)
        max_order = maximum_derivative_order(all_equations)

        spatial_derivative = SpatialDerivative(spatial_scheme, grid, max_order)

        spatial_derivatives, time_derivatives = get_derivatives(all_equations + all_formulas)

        evaluations = {}
        evaluations = create_formula_evaluations(all_formulas, evaluations)
        evaluations = create_derivative_evaluations(spatial_derivatives,evaluations, spatial_derivative)

        # We will assume that all the functions in time derivative are known at the start
        order_of_evaluations = []
        known = [d.args[0] for d in time_derivatives]
        for val in known:
            evaluated = Evaluations(val, val, None, None, val)
            evaluations[val] =  evaluated
            order_of_evaluations += [val]

        # Sort the terms in the order they should be evaluated (with respect to their dependencies).
        # First get the primitive variables that the time derivatives are applied to (e.g. u_i in Der(u_i, t))
        order_of_evaluations = sort_evaluations(order_of_evaluations, evaluations, Indexed)
        # Then sort the derivatives
        order_of_evaluations = sort_evaluations(order_of_evaluations, evaluations, Derivative)

        # Update the range of evaluations for each evaluation
        set_range_of_evaluations(order_of_evaluations, evaluations, grid)

        work_array_index = 0
        work_array_name = 'wk'
        # update the work arrays
        evaluations, work_array_index = update_work_arrays(order_of_evaluations, evaluations, work_array_name, work_array_index, grid)

        self.computations = []
        self.computations += create_formula_kernels(order_of_evaluations,evaluations, known, grid)
        derivatives = [ev for ev in order_of_evaluations if isinstance(ev, Derivative) and ev not in known]
        self.computations += create_derivative_kernels(derivatives,evaluations,\
            spatial_derivative, work_array_name, work_array_index, grid)

        # All the spatial computations are evaluated by this point. Now get the updated equations.
        updated_equations = substitute_work_arrays(order_of_evaluations,evaluations,all_equations )
        
        # The final computations of the residual (change in the RHS terms of the equations).
        # The residual equations are also named as work arrays.
        # The residual arrays are tracked for use in the evaluation of the temporal scheme.
        residual_equations = []
        residual_arrays = []
        for e in updated_equations:
            work_array = grid.work_array('%s%d' % (work_array_name, work_array_index))
            work_array_index += 1
            residual_arrays.append({e.lhs : work_array})
            residual_equations.append(Eq(work_array, e.rhs))
        evaluation_range = [tuple([0, s]) for s in grid.shape]
        self.computations.append(Kernel(residual_equations, evaluation_range, "Residual of equation", grid))
        
        # Update the residual arrays
        self.residual_arrays = residual_arrays
        
        # Vectors in the LHS of the equations, which are used in Symmetry Bc, these will be stored as
        # prognostic_classified in Temporal discretisation, use them 
        self.lhs_vectors = []
        for eq in expanded_equations:
            if len(eq) >1:
                temp = [list(e.lhs.atoms(Indexed)) for e in eq if len(list(e.lhs.atoms(Indexed))) == 1]
                self.lhs_vectors += [flatten(temp)]
            else:
                self.lhs_vectors += list(eq[0].lhs.atoms(Indexed))
                
        return

class GridBasedInitialisation(object):

    """ Initialise the equations on the grid of solution points.  """

    def __init__(self, grid, ics):
        """ Create the initialisation kernels.

        :arg grid: The numerical grid of solution points.
        :arg list ics: A list of initial condition formulas.
        :returns: None
        """

        self.computations = []
        initialisation_equation = []
        for ic in ics:
            initialisation_equation.append(parse_expr(ic, local_dict = {'grid':grid, 'Symbol': EinsteinTerm}))
        range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]

        self.computations.append(Kernel(initialisation_equation, range_of_evaluation, "Initialisation", grid))

        return


