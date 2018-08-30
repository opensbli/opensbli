
#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (c) see License file

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

"""@brief
   @authors Satya Pramod Jammy
   @contributors
   @details
"""

from sympy import flatten, Idx, Rational
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.datatypes import Int
from opensbli.schemes.spatial.scheme import Scheme, TemporalSolution
from opensbli.core.kernel import ConstantsToDeclare as CTD


class RungeKutta(Scheme):
    """ Applies a Runge-Kutta time-stepping scheme.

        :arg int order: The order of accuracy of the scheme."""
    def __init__(cls, order, constant_dt=None):
        Scheme.__init__(cls, "RungeKutta", order)
        cls.schemetype = "Temporal"
        cls.nloops = 2
        cls.stage = Idx('stage', order)
        cls.solution_coeffs = ConstantIndexed('rkold', cls.stage)
        cls.stage_coeffs = ConstantIndexed('rknew', cls.stage)
        # Update coefficient values
        cls.get_coefficients
        niter_symbol = ConstantObject('niter', integer=True)
        niter_symbol.datatype = Int()
        cls.iteration_number = Globalvariable("iter", integer=True)
        cls.iteration_number._value = None
        cls.iteration_number.datatype = Int()
        # As iteration number is used in a for loop we dont add them to constants to declare
        cls.temporal_iteration = Idx(cls.iteration_number, niter_symbol)
        CTD.add_constant(niter_symbol)
        CTD.add_constant(cls.solution_coeffs)
        CTD.add_constant(cls.stage_coeffs)
        cls.solution = {}
        if constant_dt:
            raise NotImplementedError("")
        else:
            cls.constant_time_step = True
        cls.time_step = ConstantObject("dt")
        CTD.add_constant(cls.time_step)
        return

    @property
    def get_coefficients(cls):
        """ Return the coefficients of the Runge-Kutta update equations.

        :returns: A dictionary of (update_equation, coefficients) pairs.
        :rtype: dict
        """

        if cls.order == 3:
            cls.solution_coeffs.value = [Rational(1.0, 4.0), Rational(3.0, 20), Rational(3.0, 5.0)]
            cls.stage_coeffs.value = [Rational(2, 3), Rational(5, 12), Rational(3, 5)]
        return

    def __str__(cls):
        return "%s" % (cls.__class__.__name__)

    def get_local_function(cls, list_of_components):
        from opensbli.core.opensblifunctions import TemporalDerivative
        CD_fns = []
        for c in flatten(list_of_components):
            CD_fns += list(c.atoms(TemporalDerivative))
        return CD_fns

    def discretise(cls, type_of_eq, block):
        # We need only the equations as they contain residual residual_arrays
        if type_of_eq in cls.solution.keys():
            pass
        else:
            cls.solution[type_of_eq] = TemporalSolution()
        td_fns = cls.get_local_function(type_of_eq.equations)
        # print td_fns
        if td_fns:
            # Create a Kernel for the update ()
            old_data_sets = cls.create_old_data_sets(td_fns, block)
            new_data_sets = [eq.time_advance_array for eq in td_fns]
            zipped = zip(old_data_sets, new_data_sets)

            # create a kernel for the save equations
            kernel = cls.create_start_computations(zipped, block)
            # Add Kernel to the Solution
            cls.solution[type_of_eq].start_kernels += [kernel]
            # Create the stage and solution updates
            residuals = [eq.residual for eq in flatten(type_of_eq.equations)]
            zipped = zip(old_data_sets, new_data_sets, residuals)
            kernels = cls.create_discretisation_kernel(zipped, block)
            cls.solution[type_of_eq].kernels += kernels
            # New testing
            type_of_eq.temporalsolution = TemporalSolution()
            type_of_eq.temporalsolution.kernels += kernels
            type_of_eq.temporalsolution.start_kernels += cls.solution[type_of_eq].start_kernels
        return

    def create_discretisation_kernel(cls, zipped, block):
        solution_update_kernel = Kernel(block, "Temporal solution advancement")
        stage_update_kernel = Kernel(block, "Sub stage advancement")
        # Update the range of evaluation
        solution_update_kernel.set_grid_range(block)
        stage_update_kernel.set_grid_range(block)
        # Update the solution and stages
        if cls.constant_time_step:
            old, new = cls.constant_time_step_solution(zipped)
        solution_update_kernel.add_equation(old)
        stage_update_kernel.add_equation(new)
        solution_update_kernel.update_block_datasets(block)
        stage_update_kernel.update_block_datasets(block)
        return [stage_update_kernel, solution_update_kernel]

    def constant_time_step_solution(cls, zipped):
        dt = cls.time_step
        old = []
        new = []
        for z in zipped:
            old += [OpenSBLIEq(z[0], z[0] + dt*cls.solution_coeffs*z[2], evaluate=False)]
            new += [OpenSBLIEq(z[1], z[0] + dt*cls.stage_coeffs*z[2], evaluate=False)]
        return old, new

    def create_start_computations(cls, zipped, block):
        kernel = Kernel(block)
        kernel.computation_name = "Save equations"
        kernel.set_grid_range(block)
        for z in zipped:
            kernel.add_equation(OpenSBLIEq(z[0], z[1]))
        kernel.update_block_datasets(block)
        return kernel

    def create_old_data_sets(cls, equations, block):
        old_data_sets = []
        for no, eq in enumerate(flatten(equations)):
            fn = eq.time_advance_array
            old_data_sets += [block.work_array('%s_old' % fn.base.label)]
        return old_data_sets

    def generate_inner_loop(cls, kernels):
        from opensbli.code_generation.algorithm import DoLoop
        rkloop = DoLoop(cls.stage)
        rkloop.add_components(kernels)
        return rkloop
