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

COMMENT_DELIMITER = {"OPSC": "//", "F90": "!"}
END_OF_STATEMENT_DELIMITER = {"OPSC": ";", "F90": ""}

class GenerateCode():
    ''' generates code for the given grid, solution and boundary
    '''
    tilltime_loop = []
    timeadvancecode = []
    after_timeadvance = []
    def __init__(self, grid, spatial_solution, temporal_soln, boundary, Ics, IO):

        return
class OPSCCodePrinter(CCodePrinter):
    def __init__(self, Indexed_accs):
        settings = {}
        CCodePrinter.__init__(self, settings)
        # Indexed access numbers are required in dictionary
        self.Indexed_accs = Indexed_accs
    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '(%d.0/%d.0)' %(p,q)
    def _print_Indexed(self, expr):
        # find the symbols in the indices of the expression
        syms = flatten([list(at.atoms(Symbol)) for at in expr.indices])
        # Replace the symbols in the indices with `zero'
        for x in syms:
            expr = expr.subs({x: 0})
        if self.Indexed_accs[expr.base] != None:
            out = "%s[%s(%s)]"%(self._print(expr.base.label) \
            , self.Indexed_accs[expr.base], ','.join([self._print(ind)  for ind in expr.indices]))
        else:
            out = "%s[%s]"%(self._print(expr.base.label) \
            , ','.join([self._print(ind)  for ind in expr.indices]))
        return out

def ccode(expr, Indexed_accs):
    if isinstance(expr, Eq):
        return OPSCCodePrinter(Indexed_accs).doprint(expr.lhs) \
        + ' = ' + OPSCCodePrinter(Indexed_accs).doprint(expr.rhs)
    return OPSCCodePrinter(Indexed_accs).doprint(expr)

class OPSC(object):
    '''
    '''
    # OPS Access types, used for kernel call
    ops_access = {'inputs':'OPS_READ', 'outputs':'OPS_WRITE', 'inputoutput':'OPS_RW'}
    # OPS kenel header
    ops_header = {'inputs':'const %s %s', 'outputs':'%s %s', 'inputoutput':'%s %s', 'Idx':'const int %s'}
    # Line comments
    line_comment = "//"
    # end of line delimiter
    end_of_statement = ";"
    # Commonly used brackets
    open_flower_brace = "{"; close_fower_brace = "}"
    open_brace = "("; close_brace = ")"
    def __init__(self):
        return
    def kernel_computation(self,computation, number, **assumptions):
        print "IN OPSC"
        precission = assumptions['precission']
        pprint(precission)
        out = []
        header = []
        if computation.name == None:
            name = self.get_kernel_name(number)
        name = name + OPSC.open_brace
        # process inputs
        header = ([OPSC.ops_header['inputs']%(precission,inp) for inp in computation.inputs.keys()] + \
            [OPSC.ops_header['outputs']%(precission,inp) for inp in computation.outputs.keys()] + \
                [OPSC.ops_header['inputoutput']%(precission,inp) for inp in computation.inputoutput.keys()])
        if computation.Idx:
            header += [OPSC.ops_header['Idx']%('idx') ]
        header = [name + ' , '.join(header) + OPSC.close_brace ]
        header += [OPSC.open_flower_brace]
        code = []
        ops_accs = self.get_OPS_ACCESS_number(computation)
        pprint(ops_accs)
        pprint(header)
        for eq in computation.equations:
            pprint(ccode(eq,ops_accs))
        return
    def get_OPS_ACCESS_number(self, computation):
        '''
        Returns a dictionary of OPS_ACC's
        '''
        ops_accs = {}
        allidb = list(computation.inputs.keys()) +  list(computation.outputs.keys()) + list(computation.inputoutput.keys())
        for no,inp in enumerate(allidb):
            if inp.is_grid:
                ops_accs[inp] = 'OPS_ACC%d'%no
            else:
                ops_accs[inp] = None

        return ops_accs
    def get_kernel_name(self, number):
        return 'void kernel_block_%d_computation_%d' %(number[0], number[1])


