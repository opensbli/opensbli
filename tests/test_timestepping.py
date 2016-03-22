#!/usr/bin/env python

import os
import pytest

# OpenSBLI functions
from opensbli.timestepping import *
from opensbli.grid import *

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
    
if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
