#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.problem import *
from autofd.equations import *

def test_EinsteinTerm():
    """ Check that all indices of an Einstein term are parsed correctly. """
    term = EinsteinTerm("u_j")
    assert term.indices == ["_j"]


def test_equations_to_dict():
    """ Ensure that a list of Sympy Eq objects is converted to a dictionary. """

    mass = "Eq(Der(rho,t), -Conservative(rhou_j,x_j))"
    momentum = "Eq(Der(rhou_i,t), -Conservative(rhou_i*u_j + p* KroneckerDelta(_i,_j),x_j) + Der(tau_i_j,x_j) )"
    energy = "Eq(Der(rhoE,t), -Conservative((p+rhoE)*u_j,x_j) -Der(q_i,x_i) + Der(u_i*tau_i_j ,x_j) )"
    equations = [mass, momentum, energy]

    d = equations_to_dict([e.expanded[0] for e in expanded_equations])

    assert isinstance(d, dict)

    return


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
