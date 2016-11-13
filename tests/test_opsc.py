#!/usr/bin/env python

import os
import pytest

from sympy import symbols, pi, cos

# OpenSBLI classes and functions
from opensbli.opsc import ccode

def test_ccode():
    """ Check that the OPSC code writer outputs the expected C code statement.
    For example, all ** operators should be converted to the C 'pow' function. """
    
    x, y = symbols("x y")
    expression = pi*cos(x)**y
    expected = "M_PI*pow(cos(x), y)"
    result = ccode(expression)
    assert result == expected


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
