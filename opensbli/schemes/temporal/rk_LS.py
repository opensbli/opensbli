"""@brief
   @authors David J Lusher
   @contributors
   @details
"""

from sympy import flatten, Idx, sqrt, Rational
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.datatypes import Int
from opensbli.core.grid import GridVariable
from opensbli.schemes.spatial.scheme import Scheme, TemporalSolution


class RungeKuttaLS(Scheme):
    """ Applies a low storage (2 Register, Williamson form) Runge-Kutta scheme from: Fourth-Order Kutta Schemes 2N-Storage (Carpenter, 1994).
        Optimal SSP coefficients for the 3rd order scheme are taken from "Gottlieb, Shu, Tadmor (2001): Strong Stability-Preserving High-Order Time
        Discretization Methods, SIAM Review Vol. 43, No.1, pp 89-112.

        For m stages and current time 'n' with u^0 = u^n, du^0 = 0:
        do i = 1..m:
            du^i = A[i]*du^(i-1) + dt*Residual
            u^i = u^(i-1) + B[i]*du^i
        u^n+1 = u^m."

        :arg int order: The order of accuracy of the scheme."""

    def __init__(cls, order, formulation=None, constant_dt=None):
        Scheme.__init__(cls, "RungeKutta", order)
        if order == 3:  # 3rd order schemes are 3-stage
            n_stages = order
        elif order == 4:  # 4th order scheme is 5-stage
            n_stages = order + 1
        cls.formulation = formulation
        if formulation == 'SSP':
            if order != 3:
                raise ValueError("RK-SSP is only available for the 3rd order scheme.")
            else:
                print("An SSP Runge-Kutta scheme of order %d is being used for time-stepping." % order)
        else:
            print("A Runge-Kutta scheme of order %d is being used for time-stepping." % order)
        cls.schemetype = "Temporal"
        cls.nloops = 2
        cls.stage = Idx('stage', n_stages)
        cls.solution_coeffs = ConstantIndexed('rkB', cls.stage)
        cls.stage_coeffs = ConstantIndexed('rkA', cls.stage)
        from opensbli.core.kernel import ConstantsToDeclare as CTD
        # Update coefficient values
        cls.get_coefficients
        niter_symbol = ConstantObject('niter', integer=True)
        niter_symbol.datatype = Int()
        cls.iteration_number = Globalvariable("iter", integer=True)
        cls.iteration_number._value = None
        cls.iteration_number.datatype = Int()
        # As iteration number is used in a for loop we dont add them to
        # constants to declare
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
        """ Create A (intermediate update) and B (solution advance) coefficients for the RK scheme."""
        if cls.order == 3:
            if cls.formulation == 'SSP':
                c = 0.924574
                z1 = float(sqrt(36*c**4 + 36*c**3 - 135*c**2 + 84*c - 12))
                z2 = float(2*c**2 + c - 2)
                z3 = float(12*c**4 - 18*c**3 + 18*c**2 - 11*c + 2)
                z4 = float(36*c**4 - 36*c**3 + 13*c**2 - 8*c + 4)
                z5 = float(69*c**3 - 62*c**2 + 28*c - 8)
                z6 = float(34*c**4 - 46*c**3 + 34*c**2 - 13*c + 2)
                B1 = 0.924574
                B2 = (12*c*(c-1)*(3*z2-z1) - (3*z2-z1)**2)/(144*c*(3*c-2)*(c-1)**2)
                B3 = (-24*(3*c-2)*(c-1)**2)/((3*z2-z1)**2 - 12*c*(c-1)*(3*z2-z1))
                A1 = 0.0
                A2 = (-z1*(6*c**2 - 4*c + 1) + 3*z3)/((2*c+1)*z1 - 3*(c+2)*(2*c-1)**2)
                A3 = (-z1*z4 + 108*(2*c-1)*c**5 - 3*(2*c-1)*z5)/(24*z1*c*(c-1)**4 + 72*c*z6 + 72*c**6 * (2*c-13))
                cls.solution_coeffs.value = [B1, B2, B3]
                cls.stage_coeffs.value = [A1, A2, A3]
            else:
                A1, A2, A3 = 0, Rational(-5, 9), Rational(-153, 128)
                B1, B2, B3 = Rational(1, 3), Rational(15, 16), Rational(8, 15)
                cls.solution_coeffs.value = [B1, B2, B3]
                cls.stage_coeffs.value = [A1, A2, A3]
        elif cls.order == 4:
            A1, A2, A3, A4, A5 = 0, -0.4178904745, -1.192151694643, -1.697784692471, -1.514183444257
            B1, B2, B3, B4, B5 = 0.1496590219993, 0.3792103129999, 0.8229550293869, 0.6994504559488, 0.1530572479681
            cls.solution_coeffs.value = [B1, B2, B3, B4, B5]
            cls.stage_coeffs.value = [A1, A2, A3, A4, A5]
        else:
            raise NotImplementedError("Only 3rd and 4th order RK schemes are currently implemented.")
        return

    def __str__(cls):
        return "%s" % (cls.__class__.__name__)

    def get_local_function(cls, list_of_components):
        """ Finds the time derivatives to be advanced."""

        from opensbli.core.opensblifunctions import TemporalDerivative
        CD_fns = []
        for c in flatten(list_of_components):
            CD_fns += list(c.atoms(TemporalDerivative))
        return CD_fns

    def discretise(cls, type_of_eq, block):
        """ Main discretise function for the temporal advancement."""
        # We need only the equations as they contain residual residual_arrays
        if type_of_eq in cls.solution.keys():
            pass
        else:
            cls.solution[type_of_eq] = TemporalSolution()
        td_fns = cls.get_local_function(type_of_eq.equations)
        if td_fns:
            # Create a Kernel for the update ()
            temp_data_sets = cls.create_temp_data_sets(td_fns, block)
            new_data_sets = [eq.time_advance_array for eq in td_fns]
            # Create the stage and solution updates
            residuals = [eq.residual for eq in flatten(type_of_eq.equations)]
            zipped = zip(temp_data_sets, new_data_sets, residuals)
            kernels = cls.create_discretisation_kernel(zipped, block)
            cls.solution[type_of_eq].kernels += kernels
            # New testing
            type_of_eq.temporalsolution = TemporalSolution()
            type_of_eq.temporalsolution.kernels += kernels
            type_of_eq.temporalsolution.start_kernels += cls.solution[type_of_eq].start_kernels
        return

    def create_discretisation_kernel(cls, zipped, block):
        """ Creates the kernels for the intermediate step and time update.

        :arg list zipped: List of tuples containing the intermediate, solution and residual arrays for each equation being solved.
        :arg object block: OpenSBLI SimulationBlock.
        :returns: list: List of the two discretised Kernels required for the RK scheme."""
        solution_update_kernel = Kernel(block, "Temporal solution advancement")
        # Update the range of evaluation
        solution_update_kernel.set_grid_range(block)
        # Update the solution and stages
        if cls.constant_time_step:
            solution_update = cls.constant_time_step_solution(zipped)
        solution_update_kernel.add_equation(solution_update)
        solution_update_kernel.update_block_datasets(block)
        return [solution_update_kernel]

    def constant_time_step_solution(cls, zipped):
        """ Creates the equations for the intermediate step and solution update kernels.

        :arg list zipped: List of tuples containing the intermediate, solution and residual arrays for each equation being solved.
        :returns: list solution_update: Equations for the time advancement of the solution.
        :returns: list intermediate_update: Equations for the intermediate update step."""
        dt = cls.time_step
        solution_update = []
        for z in zipped:
            solution_update += [OpenSBLIEq(z[0], cls.stage_coeffs*z[0] + dt*z[2], evaluate=False)]
            solution_update += [OpenSBLIEq(z[1], z[1] + cls.solution_coeffs*z[0], evaluate=False)]
        return solution_update

    def create_temp_data_sets(cls, equations, block):
        """Creates the arrays to store the intermediate update.

        :arg list equations: Equations to be advanced in time.
        :arg object block: OpenSBLI SimulationBlock.
        :returns: list temp_datasets: List of the temporary storage DataSets."""
        temp_data_sets = []
        for no, eq in enumerate(flatten(equations)):
            fn = eq.time_advance_array
            temp_data_sets += [block.work_array('temp_%s' % fn.base.label)]
        return temp_data_sets

    def generate_inner_loop(cls, kernels):
        from opensbli.code_generation.algorithm.algorithm import DoLoop
        rkloop = DoLoop(cls.stage)
        rkloop.add_components(kernels)
        return rkloop
