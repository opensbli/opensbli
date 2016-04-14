#!/usr/bin/env python

import os
import pytest

# OpenSBLI functions
from opensbli.array import MutableDenseNDimArray
from opensbli.grid import *


@pytest.fixture
def grid():
    return Grid(ndim=3)

@pytest.fixture
def grid_with_data():
    return Grid(ndim=3, grid_data={"number_of_points":[10, 5, 2], "delta":[0.5, 0.25, 1.0]})

@pytest.fixture
def grid_variable():
    return GridVariable('test')

def test_grid(grid, grid_with_data):
    """ Ensure that a numerical grid of solution points is set up correctly. """

    assert len(grid.shape) == 3
    assert grid.shape == tuple(symbols('nx0:3', integer=True)) # nx x ny x nz grid shape
    assert grid.indices == tuple(symbols('i0:3', integer=True)) # Grid index for each dimension
    assert len(grid.deltas) == 3

    assert grid_with_data.total_points == 10*5*2 # Grid index for each dimension
    #assert grid_with_data.deltas == [0.5, 0.25, 1.0]
    
    return


def test_work_array(grid):
    """ Ensure that a work array is set up correctly on the Grid. """
    xi = EinsteinTerm('x_i')
    indices = xi.get_array(xi.get_indexed(3)).tolist()
    assert grid.work_array("test") == IndexedBase("test")[indices]
    return

def test_grid_variable(grid, grid_variable):
    """ Ensure that a variable on a grid returns of type GridVariable. """
    
    assert grid.grid_variable('test') == grid_variable
    
    return
if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
