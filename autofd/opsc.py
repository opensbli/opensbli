#!/usr/bin/env python

#    AutoFD: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of AutoFD.

#    AutoFD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    AutoFD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with AutoFD.  If not, see <http://www.gnu.org/licenses/>.

from sympy import *
from sympy.printing.ccode import CCodePrinter
from sympy.parsing.sympy_parser import parse_expr
import re
# AutoFD functions
from .codegen_utils import END_OF_STATEMENT_DELIMITER

import logging
LOG = logging.getLogger(__name__)

class GenerateCode():
    ''' generates code for the given grid, solution and boundary
    '''
    tilltime_loop = []
    timeadvancecode = []
    after_timeadvance = []

    def __init__(self, grid, spatial_solution, temporal_soln, boundary, Ics, IO):
        
        return
class OPSCCodePrinter(CCodePrinter):
    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return str(float(p)/float(q))

def ccode(expr):
    return OPSCCodePrinter().doprint(expr)
