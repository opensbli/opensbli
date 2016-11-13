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
