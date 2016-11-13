#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

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

import os
import pytest

from sympy import Symbol, Idx, flatten

# OpenSBLI classes and functions
from opensbli.spatial import SpatialDerivative, Central
from opensbli.grid import Grid

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


def test_create_stencil(spatial_derivative, central_scheme, grid):
    """ Ensure that a central difference stencil is constructured properly on a given grid. """

    stencil = spatial_derivative.create_stencil(central_scheme, grid)
    assert len(stencil) == 2
    for dim in range(len(stencil)):
        assert len(stencil[dim]) == 5 # 5 point stencil per dimension.
        i = Symbol("i%d" % dim, integer=True)
        assert stencil[dim] == [i-2, i-1, i, i+1, i+2] # Central difference
        
    return


def test_central_scheme():
    """ Check that the points of the central scheme are computed correctly for various orders. """

    assert Central(2).points == [-1, 0, 1]
    assert Central(4).points == [-2, -1, 0, 1, 2]
    assert Central(6).points == [-3, -2, -1, 0, 1, 2, 3]
    assert Central(8).points == [-4, -3, -2, -1, 0, 1, 2, 3, 4]
        
    return


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
