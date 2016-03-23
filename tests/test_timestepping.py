#!/usr/bin/env python

import os
import pytest

from sympy import Symbol, Idx

# OpenSBLI functions
from opensbli.equations import *
from opensbli.timestepping import *
from opensbli.spatial import *
from opensbli.grid import *

@pytest.fixture
def coordinate_symbol():
    return "x"


@pytest.fixture
def mass(coordinate_symbol):
    return Equation("Eq(Der(phi,t),- c_j*Der(phi,x_j))", 2, coordinate_symbol, substitutions=[], constants=["c_j"])


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
def spatial_discretisation(mass, grid, central_scheme):
    return SpatialDiscretisation([mass.expanded], [], grid, central_scheme)


@pytest.fixture
def temporal_discretisation(rk3, grid, spatial_discretisation):
    return TemporalDiscretisation(temporal_scheme=rk3, grid=grid, const_dt=True, spatial_discretisation=spatial_discretisation)
    
    
@pytest.fixture
def rk3():
    return RungeKutta(order=3)


@pytest.fixture
def forward_euler():
    return ForwardEuler()
    
    
def test_rk3(rk3):
    """ Ensure that an RK3 time-stepping scheme is set up correctly. """

    assert rk3.order == 3
    
    # Stages of the RK scheme.
    stage = Symbol('stage', integer=True)
    old = IndexedBase('rkold')[stage]
    new = IndexedBase('rknew')[stage]
    
    coefficients = rk3.get_coefficients()
    assert coefficients[old] == [-1, 2, -1]
    assert coefficients[new] == [-1, 4, -6, 4, -1]
        
    return


def test_forward_euler(forward_euler):
    """ Ensure that the forward Euler time-stepping scheme is set up correctly. """
    assert forward_euler.order == 1
    return
    

def test_temporal_discretisation(temporal_discretisation):
    """ Ensure that the time discretisation scheme is applied correctly. """
    assert temporal_discretisation.conservative == [IndexedBase("phi")]
    
    assert len(temporal_discretisation.start_computations) == 1 # For the 'save' equations.
    assert len(temporal_discretisation.computations) == 2 # There should be 2 stages in the main body of the computation, since an RK3 scheme is being used.
    
    return
    
    
if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
