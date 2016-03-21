#!/usr/bin/env python

import os
import pytest

# OpenSBLI functions
from opensbli.problem import *

def test_expand():
    """ Ensure that an equation is expanded correctly. """

    equations = ["Eq(Der(rho,t), -c*Conservative(rhou_j,x_j))"]
    substitutions = []
    ndim = 1
    constants = ["c"]
    coordinate_symbol = "x"
    metrics = [False, False]
    formulas = ["Eq(u_i, rhou_i/rho)"]

    problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)

    expanded_equations, expanded_formulas = problem.expand()
    assert len(expanded_equations) == 1
    assert str(expanded_equations[0].expanded) == "[Derivative(rho[x0, t], t) == -c*Derivative(rhou0[x0, t], x0)]"

    assert len(expanded_formulas) == 1
    assert str(expanded_formulas[0].expanded) == "[u0[x0, t] == rhou0[x0, t]/rho[x0, t]]"

    return

if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
