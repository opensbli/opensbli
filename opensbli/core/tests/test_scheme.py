
#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (c) see License file

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

from sympy import Idx, Rational
import pytest
from opensbli.core.scheme import Scheme, CentralHalos, Central, RungeKutta
from opensbli.core.block import SimulationBlock

ndim = 3
order = 6
temporal_order = 3


@pytest.fixture
def central_class():
    return Central(order)


@pytest.fixture
def central_halos():
    return CentralHalos(order)


@pytest.fixture
def local_block():
    return SimulationBlock(ndim)


@pytest.fixture
def rk_class():
    return RungeKutta(temporal_order)


def test_central_halos(central_halos):
    # Test halos for a given order
    CH = central_halos
    assert CH.halos == [-3, 3]
    assert CH.get_halos(0) == -3
    assert CH.get_halos(1) == 3
    return


def test_scheme():
    name, order = 'test', 6
    S = Scheme(name, order)
    assert S.name == name
    assert S.order == order
    return


def test_central(central_class, local_block):
    cent = central_class
    # Test properties of the Central differencing class
    assert isinstance(cent, Scheme) is True
    assert cent.schemetype == 'Spatial'
    assert cent.points == [-3, -2, -1, 0, 1, 2, 3]
    assert isinstance(cent.halotype, CentralHalos) is True
    # Test central weight creation
    direction = 0
    weights = cent._generate_weights(direction, order, local_block)
    assert weights == [1, -6, 15, -20, 15, -6, 1]
    # Check setting of the halos on the block
    cent.set_halos(local_block)
    assert len(local_block.boundary_halos) == 3
    for halo_set in local_block.boundary_halos:
        for halos in halo_set:
            for halo in halos:
                assert isinstance(halo, CentralHalos) is True
    return


def test_runge_kutta(rk_class):
    rk = rk_class
    stage = Idx('stage', temporal_order)
    # Test properties of the RungeKutta time-stepping class
    assert isinstance(rk, Scheme) is True
    assert rk.nloops == 2
    assert rk.schemetype == 'Temporal'
    assert rk.stage == stage
    assert rk.constant_time_step is True
    # Check the 3rd order coefficients
    coeffs = rk.get_coefficients
    assert coeffs[rk.solution_coeffs] == [Rational(1.0, 4.0), Rational(3.0, 20), Rational(3.0, 5.0)]
    assert coeffs[rk.stage_coeffs] == [Rational(2, 3), Rational(5, 12), Rational(3, 5)]
    return
