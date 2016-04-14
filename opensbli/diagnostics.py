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


from sympy.tensor import Indexed
from .utils import *
from .evaluations import *
from .kernel import *

class ReductionVariable(Symbol):

    def __new__(self, var):
        self = Symbol.__xnew__(self, var)
        return self


class Reduction(object):

    def __init__(self, grid, equations, formulas, prognostic_variables, spatial_scheme, rtype, compute_every):

        self.computations = []

        all_equations = flatten(equations)
        all_formulas = flatten(formulas)

        # Get all the formulas used in the equations
        all_formulas = get_used_formulas(all_formulas, all_equations)

        spatial_derivatives, time_derivatives = get_derivatives(all_equations + all_formulas)

        evaluations = {}

        # Instance of Symbolic spatial derivative
        spatial_derivative = SymbolicDerivative(spatial_scheme, grid)

        # Create evaluations for the formulas
        evaluations = create_formula_evaluations(all_formulas, evaluations)

        # Create evaluations for the derivatives
        evaluations = create_derivative_evaluations(spatial_derivatives, evaluations, spatial_derivative)

        # Order the evaluations, such that the dependencies are evaluated first, prognostic_variables are take as known
        order_of_evaluations = []
        known = []
        for val in prognostic_variables:
            evaluations[val] = Evaluations(val, val, None, None, val)
            order_of_evaluations += [val]
            known += [val]

        # Sort formulas (Indexed Objects)
        order_of_evaluations = sort_evaluations(order_of_evaluations, evaluations, Indexed)

        # Sort derivatives
        order_of_evaluations = sort_evaluations(order_of_evaluations, evaluations, Derivative)

        # Set the range of evaluations
        set_range_of_evaluations(order_of_evaluations, evaluations, grid)

        # Work array name (This should be modified only if the name is changed in Spatial.py, to reduce the memory required for simulation
        work_array_index = 0
        work_array_name = 'wk'

        # Update the work arrays depending on the order of evaluations
        evaluations, work_array_index = update_work_arrays(order_of_evaluations, evaluations, work_array_name, work_array_index, grid)

        # Create formula computation kernels
        self.computations += create_formula_kernels(order_of_evaluations, evaluations, known, grid)
        derivatives = [ev for ev in order_of_evaluations if isinstance(ev, Derivative) and ev not in known]

        # Create Derivative computation kernels
        self.computations += create_derivative_kernels(derivatives, evaluations,
                                                       spatial_derivative, work_array_name, work_array_index, grid)

        # Get the reduction equations
        reduction_equations = self.create_reduction_equations(all_equations, rtype, grid)

        # Update the equations with the work arrays used
        update_equations = substitute_work_arrays(order_of_evaluations, evaluations, reduction_equations)

        # Create computations for the update equations
        evaluation_range = [tuple([0, s]) for s in grid.shape]
        self.computations.append(Kernel(update_equations, evaluation_range, "Reduction equations", grid))

        # Store the iteration number used to write the if statement
        self.compute_every = compute_every

        return

    def create_reduction_variables(self, equations):
        """ Create the reduction variables for the diagnostic equations.
       
        :arg equations: The diagnostic equations.
        :returns: A list with variables of type ReductionVariable.
        :rtype: list
        """
        reduction_variables = []
        for eq in equations:
            reduction_variables.append(ReductionVariable(str(eq.lhs.base)))
        return reduction_variables

    def create_reduction_equations(self, equations, rtype, grid):
        """ Create the reduction equations. The Indexed terms in the LHS of the diagnostic equations
        are changed to the type ReductionVariable and the type of reduction is applied
        Eg. summation reduction of the equation "Eq(umean[x0, x1, x2], u0[x0, x1, x2]" is written as
        umean = umean + f[x0, x1, x2]; where umean is a ReductionVariable object.
        
        :returns: A list of reduction equations
        :rtype: list
        """

        reduction_variable = self.create_reduction_variables(equations)
        reduction_equation = [None for eq in equations]

        for number, eq in enumerate(equations):
            if rtype[number] == "sum":
                reduction_equation[number] = Eq(reduction_variable[number],
                                                (eq.rhs + reduction_variable[number]))
            else:
                raise NotImplementedError("Only summation reductions are supported")

        return reduction_equation
