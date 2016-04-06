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

'''
For writing the variables that are already evaluated we are using FileIO
The various diagnostics to be supported are
1. Compute some variables and write to files (Eg. Statistics, vorticity, etc etc)
    this is sub divided into
    a. Do any computations if exist
    b. These can be after n time steps or at the end of the simulation
        i. Check point b, if after n iterations write if statement
        ii. Else at the end of simulation write it

2. Reduction of a variable i.e already evaluated or evaluation and write to a file
    this can be sub divded into
    a. Computation
    b. Reduction operation
    c. Writing to file (create a new file or append to a file)

Now the inputs to the class,
These are the equations expanded using Einstein expansion
type of Diagnostics time diagnostics, time average, reduction after n time steps,

'''

from sympy.tensor import IndexedBase, Indexed
from .utils import *
from .kernel import ReductionVariable


class Reduction():

    def __init__(self, grid, equations, formulas, prognostic_variables, spatial_scheme, rtype, compute_every):
        '''
        fname will be automatic??
        '''
        #self.type = reduction_type
        self.computations = []
        all_equations = flatten(equations)
        all_formulas = flatten(formulas)
        all_formulas = get_used_formulas(all_formulas, all_equations)
        spatial_derivatives, time_derivatives = get_derivatives(all_equations + all_formulas)
        evaluations = {}
        spatial_derivative = SymDerivative(spatial_scheme, grid)
        evaluations = create_formula_evaluations(all_formulas, evaluations)
        evaluations = create_derivative_evaluations(spatial_derivatives,evaluations, spatial_derivative)
        order_of_evaluations = []
        known = []
        for val in prognostic_variables:
            #val = self.create_indexed(val)
            evaluations[val] = Evaluations(val, val, None, None, val)
            order_of_evaluations += [val]
            known += [val]

        order_of_evaluations = sort_evaluations(order_of_evaluations, evaluations, Indexed)

        order_of_evaluations = sort_evaluations(order_of_evaluations, evaluations, Derivative)
        set_range_of_evaluations(order_of_evaluations, evaluations, grid)
        work_array_index = 0
        work_array_name = 'wk'
        # Creating kernels
        evaluations, work_array_index = update_work_arrays(order_of_evaluations, evaluations, work_array_name, work_array_index, grid)

        self.computations += create_formula_kernels(order_of_evaluations,evaluations, known, grid)
        derivatives = [ev for ev in order_of_evaluations if isinstance(ev, Derivative) and ev not in known]

        self.computations += create_derivative_kernels(derivatives,evaluations,\
            spatial_derivative, work_array_name, work_array_index, grid)
        reduction_equations = self.create_reduction_equations(all_equations, rtype, grid)
        update_equations = substitute_work_arrays(order_of_evaluations,evaluations,reduction_equations )

        # create computations for the update equations
        evaluation_range = [tuple([0, s]) for s in grid.shape]
        self.computations.append(Kernel(update_equations, evaluation_range, "Reduction equations", grid))
        self.compute_every = compute_every
        return
    def create_reduction_variables(self, equations):
        reduction_variables = []
        for eq in equations:
            reduction_variables.append(ReductionVariable(str(eq.lhs.base)))
        return reduction_variables
    def create_reduction_equations(self, equations, rtype, grid):
        """
        The LHS of the equations are Indexed Objects, this converts the LHS into reduction variable and
        applies the reduction operation
        """
        reduction_variable = self.create_reduction_variables(equations)
        reduction_equation = [None for eq in equations]
        for number, eq in enumerate(equations):
            if rtype[number] == "sum":
                reduction_equation[number] = Eq(reduction_variable[number], \
                    (eq.rhs + reduction_variable[number])/(grid.total_points))
            else:
                raise NotImplementedError("Only summation reductions are supported")
        return reduction_equation
    def create_indexed(self, var):
        '''
        This is a helper function Which will be removed after modifying spatial.py
        '''
        ind = [EinsteinTerm('x0'), EinsteinTerm('x1'), EinsteinTerm('t')]
        var = IndexedBase('%s'%var.base)
        return var[ind]
