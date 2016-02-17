#!/usr/bin/env python

import os
import pytest

from autofd.equations import *


@pytest.fixture
def mass():
    return Equation("Eq(Der(rho,t), -Conservative(rhou_i,x_j))", 2, substitutions=[], constants=[])


@pytest.fixture
def c():
    """ A scalar EinsteinTerm """
    return EinsteinTerm("c")


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


def test_has_index(u_i, tau_i_j):
    """ Check that various Einstein terms have the indices that they should have. """

    assert u_i.has_index("_i")
    assert not u_i.has_index("_j")

    assert tau_i_j.has_index("_i")
    assert tau_i_j.has_index("_j")
    assert not tau_i_j.has_index("_k")


def test_apply_index(u_i, tau_i_j):
    """ Check that an index can successfully be applied to an Einstein term. """

    new = u_i.apply_index("_j", 0)
    assert not u_i.has_index("_j")

    assert tau_i_j.has_index("_i")
    assert tau_i_j.has_index("_j")
    assert not tau_i_j.has_index("_k")


def test_equations_to_dict(mass):
    """ Ensure that a list of SymPy Eq objects is converted to a dictionary. """

    equations = [mass.expanded[0]]
    assert isinstance(equations, list)
    equations = equations_to_dict(equations)
    assert isinstance(equations, dict)


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
