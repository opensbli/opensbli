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

from sympy import flatten

# OpenSBLI classes and functions
from opensbli.problem import Problem

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

    expanded_equations = flatten(problem.get_expanded(problem.equations))
    expanded_formulas = flatten(problem.get_expanded(problem.formulas))

    assert len(expanded_equations) == 1
    assert str(expanded_equations[0].lhs) == "Derivative(rho[x0, t], t)"
    assert str(expanded_equations[0].rhs) == "-c*Derivative(rhou0[x0, t], x0)"

    assert len(expanded_formulas) == 1
    assert str(expanded_formulas[0].lhs) == "u0[x0, t]"
    assert str(expanded_formulas[0].rhs) == "rhou0[x0, t]/rho[x0, t]"

    # Test the other way of expanding equations
    assert expanded_equations == flatten(problem.get_expanded(problem.expand(equations)))

    return

if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
