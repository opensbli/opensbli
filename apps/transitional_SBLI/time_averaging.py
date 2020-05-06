from opensbli import *
from sympy import symbols, exp, pprint
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.equation_types.opensbliequations import OpenSBLIEquation as Eq
from opensbli.postprocess.post_process_eq import *
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.core.grid import GridVariable
from opensbli.code_generation.algorithm.common import *
from sympy import Function

class UserDefinedEquations(NonSimulationEquations, Discretisation, Solution):
    """User defined equations, this will not discretise the equations. No checking is performed.
    Just forms a kernel on the range and places the kernel in the algorithm place passed as an
    input to the class
    """

    def __new__(cls, **kwargs):
        ret = super(UserDefinedEquations, cls).__new__(cls)
        ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        ret._place = []
        return ret

    @property
    def algorithm_place(cls):
        return cls._place

    @algorithm_place.setter
    def algorithm_place(cls, place):
        cls._place += [place]
        return

    def spatial_discretisation(cls, block):
        """ Applies the spatial discretisation of the equations by calling the discretisation of each spatial scheme provided on the block

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """

        # Instantiate the solution class
        cls.solution = Solution()
        user_defined_kernel = Kernel(block)
        name = " ".join([p.__class__.__name__ for p in cls.algorithm_place])
        user_defined_kernel.set_computation_name("user kernel %s" % name)
        user_defined_kernel.set_grid_range(block)
        user_defined_kernel.add_equation(cls.equations)
        user_defined_kernel.update_block_datasets(block)
        cls.Kernels += [user_defined_kernel]
        return

    def apply_boundary_conditions(cls, block):
        return
    
eqns = []
niter = ConstantObject('niter')
rho, rhou0, rhou1, rhou2, rhoE = DataObject('rho'), DataObject('rhou0'), DataObject('rhou1'), DataObject('rhou2'), DataObject('rhoE')
rho_mean, rhou0_mean, rhou1_mean, rhou2_mean, rhoE_mean = DataObject('rho_mean'), DataObject('rhou0_mean'), DataObject('rhou1_mean'), DataObject('rhou2_mean'), DataObject('rhoE_mean')
rhou0u0_mean, rhou1u1_mean, rhou2u2_mean = DataObject('rhou0u0_mean'), DataObject('rhou1u1_mean'), DataObject('rhou2u2_mean')
rhou0u1_mean, rhou1u2_mean, rhou0u2_mean = DataObject('rhou0u1_mean'), DataObject('rhou1u2_mean'), DataObject('rhou0u2_mean')

stats_arrays = []
stats_arrays += [rho_mean, rhou0_mean, rhou1_mean, rhou2_mean, rhoE_mean]
stats_arrays += [rhou0u0_mean, rhou1u1_mean, rhou2u2_mean]
stats_arrays += [rhou0u1_mean, rhou1u2_mean, rhou0u2_mean]

# Sum process
def sum_statistics():
    eqns = []
    eqns += [Eq(rho_mean, rho_mean + rho)]
    eqns += [Eq(rhou0_mean, rhou0_mean + rhou0)]
    eqns += [Eq(rhou1_mean, rhou1_mean + rhou1)]
    eqns += [Eq(rhou2_mean, rhou2_mean + rhou2)]
    eqns += [Eq(rhoE_mean, rhoE_mean + rhoE)]
    rho_inv = GridVariable("rho_inv")
    eqns += [Eq(rho_inv, 1.0/rho)]
    eqns += [Eq(rhou0u0_mean, rhou0u0_mean + rhou0*rhou0*rho_inv)]
    eqns += [Eq(rhou1u1_mean, rhou1u1_mean + rhou1*rhou1*rho_inv)]
    eqns += [Eq(rhou2u2_mean, rhou2u2_mean + rhou2*rhou2*rho_inv)]
    eqns += [Eq(rhou0u1_mean, rhou0u1_mean + rhou0*rhou1*rho_inv)]
    eqns += [Eq(rhou1u2_mean, rhou1u2_mean + rhou1*rhou2*rho_inv)]
    eqns += [Eq(rhou0u2_mean, rhou0u2_mean + rhou0*rhou2*rho_inv)]
    return eqns


# Averaging process
def average_statistics():
    eqns = []
    eqns += [Eq(rho_mean, rho_mean/niter)]
    eqns += [Eq(rhou0_mean, rhou0_mean/niter)]
    eqns += [Eq(rhou1_mean, rhou1_mean/niter)]
    eqns += [Eq(rhou2_mean, rhou2_mean/niter)]
    eqns += [Eq(rhoE_mean, rhoE_mean/niter)]
    eqns += [Eq(rhou0u0_mean, rhou0u0_mean/niter)]
    eqns += [Eq(rhou1u1_mean, rhou1u1_mean/niter)]
    eqns += [Eq(rhou2u2_mean, rhou2u2_mean/niter)]
    eqns += [Eq(rhou0u1_mean, rhou0u1_mean/niter)]
    eqns += [Eq(rhou1u2_mean, rhou1u2_mean/niter)]
    eqns += [Eq(rhou0u2_mean, rhou0u2_mean/niter)]
    return eqns

def get_stats_classes():
    accumulation = UserDefinedEquations()
    accumulation.algorithm_place = InTheSimulation(frequency=False)

    sum_eq = sum_statistics()
    accumulation.add_equations(sum_eq)

    normalisation = UserDefinedEquations()
    normalisation.algorithm_place = AfterSimulationEnds()
    
    average_eq = average_statistics()
    normalisation.add_equations(average_eq)
    return [accumulation, normalisation]

def get_arrays():
    return stats_arrays
