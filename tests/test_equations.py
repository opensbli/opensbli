#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.problem import *
from autofd.equations import *

@pytest.fixture
def mass():
    return Equation("Eq(Der(rho,t), -Conservative(rhou_j,x_j))", 2, substitutions=[], constants=[])

def test_einstein_indices():
    """ Check that all indices of an Einstein term are parsed correctly. """

    term = EinsteinTerm("c")
    assert term.indices == []

    term = EinsteinTerm("u_j")
    assert term.indices == ["_j"]

    term = EinsteinTerm("tau_i_j")
    assert term.indices == ["_i", "_j"]

def test_equations_to_dict(mass):
    """ Ensure that a list of Sympy Eq objects is converted to a dictionary. """

    equations = [mass.expanded[0]]
    assert isinstance(equations, list)
    equations = equations_to_dict(equations)
    assert isinstance(equations, dict)

    return


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
