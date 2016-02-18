#!/usr/bin/env python

import os
import pytest

from autofd.codegen_utils import *
from autofd.problem import *


@pytest.fixture
def problem():
    equations = ["Eq(Der(rho,t), -c*Conservative(rhou_j,x_j))"]
    substitutions = []
    ndim = 1
    constants = ["c"]
    coordinate_symbol = ["x"]
    metrics = [False, False]
    formulas = ["Eq(u_i, rhou_i/rho)"]

    problem = Problem(equations, substitutions, ndim, constants, coordinate_symbol, metrics, formulas)
    return problem


def test_footer_code(problem):
    """ Check that the OPSC footer code is written correctly. """

    code = footer_code(problem, language="OPSC")
    assert "ops_exit();" in code[-2]
    assert "}" in code[-1]


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
