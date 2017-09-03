import os
import pytest
from opensbli.core.grid import GridVariable
from sympy import srepr, Indexed, pprint, Matrix, flatten, Equality
from sympy.core.numbers import Zero
from opensbli.core.parsing import EinsteinEquation
from opensbli.core.opensbliequations import SimulationEquations
from opensbli.core.opensbliobjects import EinsteinTerm, DataSet, DataSetBase, CoordinateObject
from opensbli.core.block import SimulationBlock
from opensbli.core.weno_opensbli import LLFWeno, Weno, EigenSystem,\
    SimpleAverage, RoeAverage, WenoHalos, LeftWenoReconstructionVariable, RightWenoReconstructionVariable, ConfigureWeno
from opensbli.physical_models.euler_eigensystem import EulerEquations
from opensbli.core.kernel import Kernel
from opensbli.core.opensblifunctions import WenoDerivative
from opensbli.physical_models.ns_physics import NSphysics

# Testing in 3D
ndim = 3
# Locations for averaging
locations = [0, 1]
# WENO scheme order and k, currently WENO3
k = 2
weno_order = 2*k-1


@pytest.fixture
def local_block():
    return SimulationBlock(ndim)


@pytest.fixture
def local_kernel():
    return Kernel(local_block(), computation_name="WENO test")


@pytest.fixture
def solution_vector():
    NS = NSphysics(ndim)
    return flatten([NS.density(), flatten(NS.momentum()), NS.total_energy()])


@pytest.fixture
def local_weno_derivatives():
    ids = solution_vector()[0].indices
    p = DataSetBase('p')[ids]
    rhoE = DataSetBase('rhoE')[ids]
    u0 = DataSetBase('u0')[ids]
    rhou0, rhou1, rhou2 = DataSetBase('rhou0')[ids], DataSetBase('rhou1')[ids], DataSetBase('rhou2')[ids]
    x0 = CoordinateObject('x0')
    x0.direction = 0
    flux = [WenoDerivative(rhou0, x0), WenoDerivative(p+rhou0*u0, x0), WenoDerivative(rhou1*u0, x0),
            WenoDerivative(rhou2*u0, x0), WenoDerivative((p+rhoE)*u0, x0)]
    return flux


@pytest.fixture
def euler_eq():
    return EulerEquations(ndim)


@pytest.fixture
def eig_system():
    ev_dict, LEV_dict, REV_dict = euler_eq().generate_eig_system()
    return EigenSystem(ev_dict, LEV_dict, REV_dict)


@pytest.fixture
def llfweno():
    return LLFWeno(eig_system().eigen_value, eig_system().left_eigen_vector,
                   eig_system().right_eigen_vector, weno_order, ndim, simple_avg())


@pytest.fixture
def simple_avg():
    return SimpleAverage(locations)


@pytest.fixture
def roe_avg():
    return RoeAverage(locations, ndim)


@pytest.fixture
def weno_main():
    return Weno(weno_order)


@pytest.fixture
def weno_config():
    return ConfigureWeno(k, 1)


@pytest.fixture
def weno_halos():
    return WenoHalos(weno_order)


@pytest.fixture
def weno_reconstruction_halos():
    return WenoHalos(weno_order, reconstruction=True)


@pytest.fixture
def local_sim_equations():
    sc1 = "**{\'scheme\':\'Weno\'}"
    mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j,%s))" % sc1
    momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s))" % sc1
    energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s))" % sc1
    eq = EinsteinEquation()
    coordinate_symbol = "x"
    base_eqns = [mass, momentum, energy]
    for i, base in enumerate(base_eqns):
        base_eqns[i] = eq.expand(base, ndim, coordinate_symbol, [], [])
    simulation_eq = SimulationEquations()
    for eqn in base_eqns:
        simulation_eq.add_equations(eqn)
    return simulation_eq


def test_LLFWeno(llfweno, solution_vector, local_weno_derivatives, local_kernel, local_sim_equations):
    # Tests Characteristic LLFCharacteristic and LLFWeno together, decouple classes later

    # Characteristic methods
    # Check pre-process equations are being stored to the provided kernel
    direction = 0
    llfweno.pre_process(direction, local_weno_derivatives, solution_vector, local_kernel)
    assert len(local_kernel.equations) is not 0
    # Check pre-process LHS are all GridVariables
    for eq in local_kernel.equations:
        assert isinstance(eq.lhs, GridVariable)
    # Check post-process LHS are all GridVariables
    for eq in local_kernel.equations:
        assert isinstance(eq.lhs, GridVariable)
    # Check WENO is a spatial scheme with direction taken from the inputs
    assert llfweno.schemetype == 'Spatial'
    assert llfweno.direction == direction
    # Check post process is creating GridVariables
    # This needs the derivatives to have reconstruction_work attribute, need to finish 
    ## Might have been from post_process error
    # for derivative in local_weno_derivatives:
    #     derivative.create_reconstruction_work_array(local_block())
    # llfweno.interpolate_reconstruction_variables(local_weno_derivatives, local_kernel)
    # llfweno.post_process(local_weno_derivatives, local_kernel)

    # Check CS Matrix creation
    req_points = sorted(list(set(llfweno.reconstruction_classes[0].func_points + llfweno.reconstruction_classes[1].func_points)))
    stencil, mat = llfweno.solution_vector_to_characteristic(solution_vector, direction, 'test')
    # Check the shape of CS matrix and their datatypes
    assert mat.shape == (ndim+2, len(req_points))
    for item in flatten(mat):
        assert isinstance(item, GridVariable) is True
    # Check only the required points are being used in the CS equations
    for item in flatten(stencil):
        for dsets in item.atoms(DataSet):
            assert dsets.indices[direction] in req_points
    # Check CF Matrix creation
    stencil, mat = llfweno.flux_vector_to_characteristic(local_weno_derivatives, direction, 'test')
    # Check the shape of CF matrix and their datatypes
    assert mat.shape == (ndim+2, len(req_points))
    for item in flatten(mat):
        assert isinstance(item, GridVariable) is True
    # Check only the required points are being used in the CS equations
    for item in flatten(stencil):
        for dsets in item.atoms(DataSet):
            assert dsets.indices[direction] in req_points

    # LLFCharacteristic methods
    # Check wave speed grid variables are being created, there should be ndim+2 of them
    mat, eqns = llfweno.create_max_characteristic_wave_speed(local_kernel.equations, direction)
    mat = [i for i in flatten(mat) if isinstance(i, GridVariable)]
    assert len(mat) == ndim + 2
    # Check flux splitting is being used
    assert llfweno.flux_split is True

    # LLFWeno methods
    # Check WENO grouping is working correctly
    input_eqns = flatten(local_sim_equations.equations[:])
    grouped_eqns = llfweno.group_by_direction(input_eqns)
    # Check all directions are present
    assert grouped_eqns.keys() == [i for i in range(ndim)]
    # Check the grouping only has derivatives for that direction and that they are WenoDerivatives of correct number
    for dire in [i for i in range(ndim)]:
        assert len(grouped_eqns[dire]) == ndim + 2
        for derivative in grouped_eqns[dire]:
            assert isinstance(derivative, WenoDerivative) is True
            assert isinstance(derivative.args[-1], CoordinateObject) is True
            assert derivative.args[-1].direction == dire
    # Check WENO reconstruction halos are [-1, 1]
    assert llfweno.reconstruction_halotype(weno_order).halos == [-1, 1]

    # Need to fix this
    # llfweno.discretise(local_sim_equations, local_block())
    return


def test_halos(weno_halos, weno_reconstruction_halos):
    # Check halo values
    assert weno_halos.halos == [-k, k+1]
    assert weno_reconstruction_halos.halos == [-1, 1]
    return


def test_WenoConfig(weno_config):
    # Check there are the correct number of weights
    assert len(weno_config.opt_weights) == k
    assert len(weno_config.c_rj) == k**2
    # Check symbolic dictionary is being created correctly
    points, sym_dict = weno_config.generate_symbolic_function_points
    for index, fn in sym_dict.iteritems():
        assert fn.args[-1] == index
        assert fn in points
        assert isinstance(fn, Indexed) is True
    # Check optimal weights sum to 1 as required
    assert sum([d for d in weno_config.opt_weights.values()]) == 1
    return


def test_Weno(weno_main):
    # Check order
    assert weno_main.order == 3
    assert weno_main.k == 2
    # Check halo type
    assert isinstance(weno_main.halotype, WenoHalos) is True
    # Check right and left reconstructions are present in the [right, left] order
    assert isinstance(weno_main.reconstruction_classes[0], RightWenoReconstructionVariable)
    assert isinstance(weno_main.reconstruction_classes[1], LeftWenoReconstructionVariable)
    # Checking attributes of the WENO reconstruction variables
    right = weno_main.reconstruction_classes[0]
    left = weno_main.reconstruction_classes[1]
    assert len(right.alpha_evaluated) == k
    assert len(right.alpha_symbols) == k
    assert len(right.omega_evaluated) == k
    assert len(right.omega_symbols) == k
    assert len(right.smoothness_indicators) == k
    assert len(right.smoothness_symbols) == k
    assert len(right.stencil_points) == weno_order
    assert len(right.func_points) == weno_order
    assert len(left.alpha_evaluated) == k
    assert len(left.alpha_symbols) == k
    assert len(left.omega_evaluated) == k
    assert len(left.omega_symbols) == k
    assert len(left.smoothness_indicators) == k
    assert len(left.smoothness_symbols) == k
    assert len(left.stencil_points) == weno_order
    assert len(left.func_points) == weno_order
    return


def test_RoeAverage(roe_avg):
    required = [EinsteinTerm('rho'), EinsteinTerm('a')] + [EinsteinTerm('u%d' % i) for i in range(ndim)]
    eqns = roe_avg.average(required, 1, 'AVG')
    # Check equations are being created with only [0,1] locations
    locs = [DataSetBase('u')[[0, locations[0], 0]].indices, DataSetBase('u')[[0, locations[1], 0]].indices]
    for eq in eqns:
        assert isinstance(eq, Equality)
        for dset in eq.atoms(DataSet):
            assert dset.indices in locs
        # Check LHS are all GridVariable
        assert isinstance(eq.lhs, GridVariable) is True
    # Check the output of get_dsets function
    dset1, dset2 = roe_avg.get_dsets(DataSetBase('u'))
    assert dset1.indices == locs[0]
    assert dset2.indices == locs[1]
    return


def test_SimpleAverage(simple_avg):
    required = [EinsteinTerm('rho'), EinsteinTerm('a')] + [EinsteinTerm('u%d' % i) for i in range(ndim)]
    name = 'test'
    eqns = simple_avg.average(required, 1, name)
    # Check equations are being created with only [0,1] locations
    locs = [DataSetBase('u')[[0, locations[0], 0]].indices, DataSetBase('u')[[0, locations[1], 0]].indices]
    for eq in eqns:
        assert isinstance(eq, Equality)
        for dset in eq.atoms(DataSet):
            assert dset.indices in locs
        # Check LHS are all GridVariable
        assert isinstance(eq.lhs, GridVariable) is True
    return


def test_EulerEigensystems(euler_eq):
    # Check eigensystems
    keys = [i for i in range(ndim)]
    ev_dict, LEV_dict, REV_dict = euler_eq.generate_eig_system()
    # Check all directions are present in the eigensystem dictionaries
    assert keys == list(ev_dict.keys())
    collected_terms = set()
    # Check the type and shape of the matrices
    for item in [ev_dict, LEV_dict, REV_dict]:
        for i in range(ndim):
            assert type(item[i]) == Matrix
            assert item[i].shape == (ndim+2, ndim+2)
            for variable in flatten(item[i]):
                for ET in variable.atoms():
                    if type(ET) == EinsteinTerm:
                        collected_terms.add(ET)
    # Check only a, rho and velocity components are present in the eigensysytems
    # NOTE: this will change when metrics are included in the eigensystem
    expected_terms = [EinsteinTerm('a'), EinsteinTerm('rho')] + [EinsteinTerm('u%d' % i) for i in range(ndim)]
    assert set(expected_terms) == collected_terms
    return


def test_EigenSystem(eig_system):
    # Check the correct terms are in the EigenSystem for 3D
    assert eig_system.get_symbols_in_ev(0) == set([EinsteinTerm('a'), EinsteinTerm('u0')])
    assert eig_system.get_symbols_in_ev(1) == set([EinsteinTerm('a'), EinsteinTerm('u1')])
    assert eig_system.get_symbols_in_ev(2) == set([EinsteinTerm('a'), EinsteinTerm('u2')])
    # Check the correct terms are in the LEV matrix
    full_set = set([EinsteinTerm('a'), EinsteinTerm('u0'), EinsteinTerm('u1'), EinsteinTerm('u2'), EinsteinTerm('rho')])
    assert eig_system.get_symbols_in_LEV(0) == full_set
    assert eig_system.get_symbols_in_LEV(1) == full_set
    assert eig_system.get_symbols_in_LEV(2) == full_set
    # Check the correct terms are in the REV matrix
    assert eig_system.get_symbols_in_REV(0) == full_set
    assert eig_system.get_symbols_in_REV(1) == full_set
    assert eig_system.get_symbols_in_REV(2) == full_set
    # Check GridVariables and zero elements are being generated correctly
    grid_EV = eig_system.generate_grid_variable_ev(1, 'test')
    grid_LEV = eig_system.generate_grid_variable_LEV(1, 'test')
    grid_REV = eig_system.generate_grid_variable_REV(1, 'test')
    matrix_elements = flatten(grid_EV) + flatten(grid_LEV) + flatten(grid_REV)
    for element in matrix_elements:
        assert (isinstance(element, GridVariable) or isinstance(element, Zero)) is True
    # Check EinstinTerm to GridVariable conversion
    input_matrix = eig_system.left_eigen_vector[0]
    output_matrix = eig_system.convert_matrix_to_grid_variable(input_matrix, 'test')
    # Check output shape is the same as the input
    assert output_matrix.shape == input_matrix.shape
    # No EinsteinTerms should remain
    for element in flatten(output_matrix):
        for item in element.atoms():
            assert type(item) is not EinsteinTerm
    # Check equations are being created correctly
    equations = eig_system.generate_equations_from_matrices(input_matrix, output_matrix)
    for eqn in flatten(equations):
        assert (isinstance(eqn, Equality) or isinstance(eqn, Zero)) is True
    # Check DataSet creation from an expression has the correct indices
    expr_matrix = eig_system.eigen_value[1]
    output = eig_system.convert_symbolic_to_dataset(expr_matrix, 3, 0)
    assert output.shape == expr_matrix.shape
    terms = [element for element in flatten(output) if type(element) is not Zero]
    expected_indices = [3, 0, 0]
    for expr in terms:
        for dset in expr.atoms(DataSet):
            assert len(dset.indices) == ndim + 1
            for i in range(ndim):
                assert dset.indices[i] == expected_indices[i]
    return


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
