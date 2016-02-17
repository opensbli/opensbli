#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.problem import *
from autofd.equations import *

@pytest.fixture
def mass():
    return Equation("Eq(Der(rho,t), -Conservative(rhou_j,x_j))", 2, substitutions=[], constants=[])


def test_indices():
    """ Check that all indices of an Einstein term are parsed correctly. """

    assert EinsteinTerm("c").indices == []  # Scalar
    assert EinsteinTerm("u_j").indices == ["_j"]  # Vector
    assert EinsteinTerm("tau_i_j").indices == ["_i", "_j"]  # Tensor


def test_name():
    """ Check that the name of an Einstein term is correct. """

    assert EinsteinTerm("c").name == "c"  # Scalar
    assert EinsteinTerm("u_j").name == "u_j"  # Vector
    assert EinsteinTerm("tau_i_j").name == "tau_i_j"  # Tensor


def test_equations_to_dict(mass):
    """ Ensure that a list of SymPy Eq objects is converted to a dictionary. """

    equations = [mass.expanded[0]]
    assert isinstance(equations, list)
    equations = equations_to_dict(equations)
    assert isinstance(equations, dict)


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
