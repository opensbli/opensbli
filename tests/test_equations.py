#!/usr/bin/env python

import os
import pytest

from autofd.equations import *


@pytest.fixture
def mass():
    return Equation("Eq(Der(rho,t), -Conservative(rhou_i,x_j))", 3, substitutions=[], constants=[])


@pytest.fixture
def c():
    """ A scalar EinsteinTerm """
    return EinsteinTerm("c")


@pytest.fixture
def c_constant():
    """ A constant scalar EinsteinTerm """
    c = EinsteinTerm("c")
    c.constant = True
    return c


@pytest.fixture
def u_i():
    """ A vector EinsteinTerm """
    return EinsteinTerm("u_i")


@pytest.fixture
def tau_i_j():
    """ A tensor EinsteinTerm """
    return EinsteinTerm("tau_i_j")


def test_indices(c, u_i, tau_i_j):
    """ Check that all indices of an Einstein term are parsed correctly. """

    assert c.indices == []  # Scalar
    assert u_i.indices == ["_i"]  # Vector
    assert tau_i_j.indices == ["_i", "_j"]  # Tensor


def test_name(c, u_i, tau_i_j):
    """ Check that the name of an Einstein term is correct. """

    assert c.name == "c"  # Scalar
    assert u_i.name == "u_i"  # Vector
    assert tau_i_j.name == "tau_i_j"  # Tensor


def test_has_index(c, u_i, tau_i_j):
    """ Check that various Einstein terms have the indices that they should have. """

    assert not c.has_index("_l")

    assert u_i.has_index("_i")
    assert not u_i.has_index("_j")

    assert tau_i_j.has_index("_i")
    assert tau_i_j.has_index("_j")
    assert not tau_i_j.has_index("_k")


def test_expand_index(c_constant, u_i, tau_i_j):
    """ Check that an index can successfully be applied to an Einstein term. """

    new = c_constant.expand_index("_i", 0)
    assert new == c_constant  # This should stay the same since the Einstein term is constant.

    new = u_i.expand_index("_i", 0)
    assert not new.has_index("_i")

    new = tau_i_j.expand_index("_i", 0)
    assert new.has_index("_j")
    assert not new.has_index("_i")


def test_equations_to_dict(mass):
    """ Ensure that a list of SymPy Eq objects is converted to a dictionary. """

    equations = [mass.expanded[0]]
    assert isinstance(equations, list)
    equations = equations_to_dict(equations)
    assert isinstance(equations, dict)


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
