import os
import pytest
from opensbli.core.grid import GridVariable
from sympy import Idx, srepr, Indexed, pprint, Matrix, flatten, Symbol, Equality
from sympy.core.numbers import Zero
from opensbli.core.opensbliobjects import EinsteinTerm, DataSet
from opensbli.core.block import SimulationBlock
from opensbli.core.weno_opensbli import LLFCharacteristic, LLFWeno, EigenSystem, Characteristic, SimpleAverage, RoeAverage
from opensbli.physical_models.euler_eigensystem import EulerEquations

ndim = 3


@pytest.fixture
def local_block():
    return SimulationBlock(ndim, blocknumber=0)


@pytest.fixture
def euler_eq():
    return EulerEquations(ndim)


@pytest.fixture
def eig_system():
    ev_dict, LEV_dict, REV_dict = euler_eq().generate_eig_system()
    return EigenSystem(ev_dict, LEV_dict, REV_dict)


@pytest.fixture
def characteristic():
    ev_dict, LEV_dict, REV_dict = euler_eq().generate_eig_system()
    return Characteristic(ev_dict, LEV_dict, REV_dict)


@pytest.fixture
def simple_avg():
    return SimpleAverage([0, 1])


@pytest.fixture
def roe_avg():
    return RoeAverage([0, 1], ndim)


def test_RoeAverage(roe_avg):
    return


def test_SimpleAverage(simple_avg):
    return


def test_Characteristic(characteristic):
    return


def test_LLFCharacteristic():
    return


def test_LLFWeno():
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
