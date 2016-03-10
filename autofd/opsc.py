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
from .equations import equations_to_dict

import logging
LOG = logging.getLogger(__name__)

class OPSC(object):
    def __init__(self):
        return

    def generate(self, equations, system):
        all_calls = []
        all_kernels = []
        if isinstance(equations, dict):
            for key, value in equations.iteritems():
                if isinstance(value, list):
                    call, kernel = self.get_kernel(value, system)
                    all_calls = all_calls + [call]
                    all_kernels = all_kernels + [kernel]
                else:
                    call, kernel = self.get_kernel([value], system)
                    all_calls = all_calls + [call]
                    all_kernels = all_kernels + [kernel]
        elif isinstance(equations, list):
            call, kernel = self.get_kernel(equations, system)
            all_calls = all_calls + [call]
            all_kernels = all_kernels + [kernel]
        else:
            call, kernel = self.get_kernel([equations], system)
            all_calls = all_calls + [call]
            all_kernels = all_kernels + [kernel]

        return all_calls, all_kernels

    def get_kernel(self, equations, system):
        lhs = flatten(list(list(eq.lhs.atoms(Indexed)) for eq in equations))
        rhs = flatten(list(list(eq.rhs.atoms(Indexed)) for eq in equations))
        all_indexed = list(set(lhs+rhs))

        lhs_base_labels = set(list(i.base.label for i in lhs))
        rhs_base_labels = set(list(i.base.label for i in rhs))
        inouts = list(lhs_base_labels.intersection(rhs_base_labels))
        ins = list(rhs_base_labels.difference(inouts))
        outs = list(lhs_base_labels.difference(inouts))
        all_base = ins+outs+inouts

        symbols = flatten(list(list(eq.atoms(Symbol)) for eq in equations))
        Idxs = list(set(list(i.indices for i in (lhs+rhs))))
        Idxs = flatten(list(i) for i in Idxs)
        index_labels = list(set(list(i.label for i in Idxs)))
        for i in Idxs:
            index_labels = index_labels + [i.lower, i.upper]
        symbols = list(set(symbols).difference(set(index_labels)).difference(set(system.constants)).difference(set(all_base)))
        symbols = list(set(symbols))
        out = []
        symbol_declaration = []
        equations_dict = equations_to_dict(equations)
        for s in symbols:
            # Ignore the IDX terms - these are special terms representing each grid point's index, and will be replaced with idx[0], idx[1] and idx[2] later to ensure that the idx array is used in the kernel.
            if str(s).startswith("IDX"):
                continue
            symbol_declaration = symbol_declaration + ['double %s;' % s]
            if equations_dict.get(s):
                pass
            else:
                raise ValueError("I don't know the formula for %s" % s)

        # Grid range
        lower = []
        upper = []
        for dim in range(system.ndim):
            lower = lower + [system.blockdims[dim].lower - system.halos[dim]]
            upper = upper + [system.blockdims[dim].upper + system.halos[dim]+1]

        # Create the C code representing each equation.
        for e in equations:
            code = ccode(e)
            code = code.replace('==', '=')
            code = code.replace('IDX0', 'idx[0]').replace('IDX1', 'idx[1]').replace('IDX2', 'idx[2]')
            code += END_OF_STATEMENT_DELIMITER['OPSC']
            out = out + [code]

        kernel_calls = []
        kernel_header = []
        kernel = []
        kernel_name = system.kernel_name % system.kernel_index
        system.kernel_index = system.kernel_index + 1
        kernel_call = 'ops_par_loop(%s, \"%s\", %s[%s], %d, %%s' % (kernel_name, kernel_name, system.block_name, system.block, system.ndim)
        head = 'void %s(' % kernel_name

        # First handle the array of all the grid point indices ("idx")
        kernel_header = kernel_header + ["const int *idx"]
        kernel_calls = kernel_calls + ["ops_arg_idx()"]

        if all_base:
            for ind, v in enumerate(all_base):
                if v in ins:
                    opstype = 'OPS_READ'
                    headty = 'const double *%s'
                elif v in outs:
                    opstype = 'OPS_WRITE'
                    headty = 'double *%s'
                elif v in inouts:
                    opstype = 'OPS_RW'
                    headty = 'double *%s'
                else:
                    raise ValueError("Don't know what the base is %s" % v)
                varib = flatten(list(v1 for v1 in all_indexed if v1.base.label == v))
                varib = list(set(varib))
                variabind = flatten(list(v1.indices) for v1 in varib)
                variabind = list(set(variabind))

                if all(va.upper == system.block.upper for va in variabind):
                    indexes = list(va for va in variabind)
                    for dim in range(system.ndim):
                        indexes = list(str(te).replace('x%d' % dim, '0') for te in indexes)
                    indexes = list(parse_expr(v) for v in indexes)
                    for inde in range(len(indexes)):
                        for ou in range(len(out)):
                            if isinstance(indexes[inde], tuple):
                                new = '%s[OPS_ACC%d%s]' % (v, ind, indexes[inde])
                            else:
                                new = '%s[OPS_ACC%d(%s)]' % (v, ind, indexes[inde])
                            old = ('%s\[%s\]' % (v, variabind[inde])).replace('+', '\+')
                            out[ou] = re.sub(r"\b(%s)" % old, new, out[ou])

                    # Get the stencil name to be written
                    indexes = indexes + [parse_expr(', '.join(list(str(0) for dim in range(system.ndim))))]
                    indexes = list(set(indexes))
                    if system.ndim > 1:
                        for dim in range(system.ndim):
                            indexes = sorted(indexes, key=lambda indexes: indexes[dim])
                        temp = flatten(list(list(t) for t in indexes))
                    else:
                        indexes = [sorted(indexes)]
                        temp = flatten(list(t) for t in indexes)

                    stencil = ','.join(list(str(t) for t in temp))
                    if system.stencils.get(stencil):
                        stencil_name = system.stencils.get(stencil)
                    else:
                        stencil_name = system.stencil_name % system.stencil_index
                        system.stencils[stencil] = stencil_name
                        system.stencil_index = system.stencil_index + 1

                    # Update range over which the loop is to be iterated
                    if len(indexes) == 1:
                        for dim in range(system.ndim):
                            lower[dim] = lower[dim] - indexes[0][dim]
                            upper[dim] = upper[dim] - indexes[0][dim]
                    else:
                        for dim in range(system.ndim):
                            lower[dim] = lower[dim] - indexes[0][dim]
                            upper[dim] = upper[dim] - indexes[-1][dim]
                    datatype = 'double'
                    arg_call = '%%s(%%s[%s], 1, %%s, \"%%s\", %%s)' % system.block
                    call = arg_call % ('ops_arg_dat', v, stencil_name, datatype, opstype)
                    kernel_calls = kernel_calls + [call]
                    kernel_header = kernel_header + [headty % v]
                else:
                    indexes = list(va for va in variabind)
                    indexes = list(str(te).replace(str(te), '0') for te in indexes)
                    indexes = list(parse_expr(v) for v in indexes)

                    for inde in range(len(indexes)):
                        for ou in range(len(out)):
                            temp = [indexes[inde]]
                            temp = list(str(te) for te in temp)
                            new = '%s[%s]' % (v, ','.join(temp))
                            old = str(varib[inde])
                            out[ou] = out[ou].replace(old, new)
                    datatype = 'double'
                    arg_call = '%%s(&%%s[%s], 1, \"%%s\", %%s)' % variabind[0]
                    call = arg_call % ('ops_arg_gbl', v, datatype, opstype)
                    kernel_calls = kernel_calls + [call]
                    kernel_header = kernel_header + [headty % v]
            iter_range = []
            for dim in range(system.ndim):
                iter_range = iter_range + [str(lower[dim])] + [str(upper[dim])]
            iter_range = ','.join(iter_range)
            kernel_calls.insert(0, kernel_call % 'iter_range%d' % system.iterrange)
            for indno in range(len(kernel_calls)-1):
                kernel_calls[indno] = kernel_calls[indno] + ','
            kernel_calls[-1] = kernel_calls[-1] + ');'
            kernel_calls = ['int iter_range%d[] = {%s};\n' % (system.iterrange, iter_range)] + kernel_calls
            system.iterrange = system.iterrange + 1
            kernel_header = head + ', '.join(kernel_header) + '){'
            kernel = [kernel_header] + symbol_declaration + out + ['}']
        else:
            LOG.debug(all_base)
            pass
        return kernel_calls, kernel

    

class OPSCCodePrinter(CCodePrinter):
    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return str(float(p)/float(q))

def ccode(expr):
    return OPSCCodePrinter().doprint(expr)
