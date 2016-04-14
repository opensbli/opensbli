#!/usr/bin/env python

import os
import pytest

# OpenSBLI functions
from opensbli.spatial import *
from opensbli.grid import *
from opensbli.problem import *
from opensbli.utils import *

@pytest.fixture
def grid():
    return Grid(ndim=2)


@pytest.fixture
def central_scheme():
    return Central(order=4)
    
    
@pytest.fixture
def max_order():
    return 1
    
    
@pytest.fixture
def spatial_derivative(grid, central_scheme, max_order):
    return SpatialDerivative(central_scheme, grid, max_order)


@pytest.fixture
def spatial_discretisation(navier_stokes_problem, grid, central_scheme):
    expanded_equations = navier_stokes_problem.get_expanded(navier_stokes_problem.equations)
    expanded_formulas = navier_stokes_problem.get_expanded(navier_stokes_problem.formulas)
    return SpatialDiscretisation(expanded_equations, expanded_formulas, grid, central_scheme)


@pytest.fixture
def coordinate_symbol():
    return "x"


@pytest.fixture
def mass(coordinate_symbol):
    return Equation("Eq(Der(rho,t), -Conservative(rhou_j,x_j))", 2, coordinate_symbol, substitutions=[], constants=[])



def test_create_stencil(spatial_derivative, central_scheme, grid):
    """ Ensure that a central difference stencil is constructured properly on a given grid. """

    stencil = spatial_derivative.create_stencil(central_scheme, grid)
    assert len(stencil) == 2
    for dim in range(len(stencil)):
        assert len(stencil[dim]) == 5 # 5 point stencil per dimension.
        i = Symbol("i%d" % dim, integer=True)
        assert stencil[dim] == [i-2, i-1, i, i+1, i+2] # Central difference
        
    return


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
