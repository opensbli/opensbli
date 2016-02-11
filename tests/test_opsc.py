#!/usr/bin/env python

import os
import pytest

# AutoFD functions
from autofd.opsc import *

def test_ccode():
    """ Check that the OPSC code writer outputs the expected C code statement.
    For example, all ** operators should be converted to the C 'pow' function. """
    
    x, y = symbols("x y")
    expression = pi*cos(x)**y
    expected = "M_PI*pow(cos(x), y)"
    result = ccode(expression)
    assert result == expected


def test_print_rational():
    """ Check that the OPSC code writer converts all rationals to their final floating point representation. """
    x = symbols("x")
    expression = Rational("2.0/3.0")
    expected = "0.666666666667"
    opsc = OPSCCodePrinter()
    result = opsc._print_Rational(expression)
    assert result == expected


if __name__ == '__main__':
    pytest.main(os.path.abspath(__file__))
