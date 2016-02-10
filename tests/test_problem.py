#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.problem import *

def test_expand():
    """ Ensure that an equation is expanded correctly. """

    equations = ["Eq(Der(rho,t),- c*conser(rhou_j,x_j))"]
    substitutions = []
    ndim = 1
    constants = ["c"]
    coordinate_symbol = ["x"]
    metrics = [False, False]
    formulas = ["Eq(u_i, rhou_i/rho)"]

    problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)

    expanded_equations, expanded_formulas = problem.expand()
    assert len(expanded_equations) == 1
    assert str(expanded_equations[0].expanded) == "[Derivative(rho, t) == -c*Derivative(rhou0, x0)]"
    assert str(expanded_equations[0].constants) == "[c]"

    assert len(expanded_formulas) == 1
    assert str(expanded_formulas[0].expanded) == "[u0 == rhou0/rho]"
    assert str(expanded_formulas[0].constants) == "[c]"

    return

if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
