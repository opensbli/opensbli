#!/usr/bin/env python

import os
import pytest

# OpenSBLI functions
from opensbli.array import MutableDenseNDimArray
from opensbli.grid import *


@pytest.fixture
def grid():
    return Grid(ndim=3)


def test_grid(grid):
    """ Ensure that a numerical grid of solution points is set up correctly. """

    assert len(grid.shape) == 3
    assert grid.shape == tuple(symbols('nx0:3', integer=True)) # nx x ny x nz grid shape
    assert grid.indices == tuple(symbols('i0:3', integer=True)) # Grid index for each dimension
    assert len(grid.deltas) == 3
    
    return


def test_work_array(grid):
    """ Ensure that a work array is set up correctly on the Grid. """
    indices = tuple(symbols('i0:3', integer=True))
    assert grid.work_array("test") == IndexedBase("test")[indices]
    return
    
    
if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
