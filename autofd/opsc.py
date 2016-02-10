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
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_application)
transformations = standard_transformations + (implicit_application,)
import re

# AutoFD functions
from .codegen_utils import END_OF_STATEMENT_DELIMITER
from .equations import equations_to_dict

import logging
LOG = logging.getLogger(__name__)

def OPSC_write_kernel(eqs, inp):
    def get_kernel(evals):
        lhs = flatten(list(list(eq.lhs.atoms(Indexed)) for eq in evals))
        rhs = flatten(list(list(eq.rhs.atoms(Indexed)) for eq in evals))
        all_indexed = list(set(lhs+rhs))

        lhs_base_labels = set(list(i.base.label for i in lhs))
        rhs_base_labels = set(list(i.base.label for i in rhs))
        inouts = list(lhs_base_labels.intersection(rhs_base_labels))
        ins = list(rhs_base_labels.difference(inouts))
        outs = list(lhs_base_labels.difference(inouts))
        all_base = ins+outs+inouts

        symbols = flatten(list(list(eq.atoms(Symbol)) for eq in evals))
        Idxs = list(set(list(i.indices for i in (lhs+rhs))))
        Idxs = flatten(list(i) for i in Idxs)
        index_labels = list(set(list(i.label for i in Idxs)))
        for i in Idxs:
            index_labels = index_labels + [i.lower, i.upper]
        symbols = list(set(symbols).difference(set(index_labels)).difference(set(inp.constants)).difference(set(all_base)))
        symbols = list(set(symbols))
        out = []
        symbol_declaration = []
        evdict = equations_to_dict(evals)
        for s in symbols:
            symbol_declaration = symbol_declaration + ['double %s;' % s]
            if evdict.get(s):
                pass
            else:
                raise ValueError("I don't know the formula for %s" % s)

        # Grid range
        lower = []
        upper = []
        for dim in range(inp.ndim):
            lower = lower + [inp.blockdims[dim].lower - inp.halos[dim]]
            upper = upper + [inp.blockdims[dim].upper + inp.halos[dim]+1]

        for e in evals:
            code = ccode(e)
            code = code.replace('==', '=') + END_OF_STATEMENT_DELIMITER['OPSC']
            out = out + [code]

        kernel_calls = []
        kernel_header = []
        kernel = []
        kernel_name = inp.kernel_name % inp.kernel_index
        inp.kernel_index = inp.kernel_index + 1
        kernel_call = 'ops_par_loop(%s, \"%s\", %s[%s], %d, %%s' % (kernel_name, kernel_name, inp.block_name, inp.block, inp.ndim)
        head = 'void %s(' % kernel_name
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

                if all(va.upper == inp.block.upper for va in variabind):
                    indexes = list(va for va in variabind)
                    for dim in range(inp.ndim):
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
                    indexes = indexes + [parse_expr(', '.join(list(str(0) for dim in range(inp.ndim))))]
                    indexes = list(set(indexes))
                    if inp.ndim > 1:
                        for dim in range(inp.ndim):
                            indexes = sorted(indexes, key=lambda indexes: indexes[dim])
                        temp = flatten(list(list(t) for t in indexes))
                    else:
                        indexes = [sorted(indexes)]
                        temp = flatten(list(t) for t in indexes)

                    stencil = ','.join(list(str(t) for t in temp))
                    if inp.stencils.get(stencil):
                        stencil_name = inp.stencils.get(stencil)
                    else:
                        stencil_name = inp.stencil_name % inp.stencil_index
                        inp.stencils[stencil] = stencil_name
                        inp.stencil_index = inp.stencil_index + 1

                    # Update range on which the loop to be iterated
                    if len(indexes) == 1:
                        for dim in range(inp.ndim):
                            lower[dim] = lower[dim] - indexes[0][dim]
                            upper[dim] = upper[dim] - indexes[0][dim]
                    else:
                        for dim in range(inp.ndim):
                            lower[dim] = lower[dim] - indexes[0][dim]
                            upper[dim] = upper[dim] - indexes[-1][dim]
                    datatype = 'double'
                    arg_call = '%%s(%%s[%s], 1, %%s, \"%%s\", %%s)' % inp.block
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
            for dim in range(inp.ndim):
                iter_range = iter_range + [str(lower[dim])] + [str(upper[dim])]
            iter_range = ','.join(iter_range)
            kernel_calls.insert(0, kernel_call % 'iter_range%d' % inp.iterrange)
            for indno in range(len(kernel_calls)-1):
                kernel_calls[indno] = kernel_calls[indno] + ','
            kernel_calls[-1] = kernel_calls[-1] + ');'
            kernel_calls = ['int iter_range%d[] = {%s};\n' % (inp.iterrange, iter_range)] + kernel_calls
            inp.iterrange = inp.iterrange + 1
            kernel_header = head + ', '.join(kernel_header) + '){'
            kernel = [kernel_header] + symbol_declaration + out + ['}']
        else:
            LOG.debug(all_base)
            pass
        return kernel_calls, kernel

    allcalls = []
    allkernels = []
    if isinstance(eqs, dict):
        for key, value in eqs.iteritems():
            if isinstance(value, list):
                call, comp = get_kernel(value)
                allcalls = allcalls + [call]
                allkernels = allkernels + [comp]
            else:
                call, comp = get_kernel([value])
                allcalls = allcalls + [call]
                allkernels = allkernels + [comp]
    elif isinstance(eqs, list):
        call, comp = get_kernel(eqs)
        allcalls = allcalls + [call]
        allkernels = allkernels + [comp]
    else:
        call, comp = get_kernel([eqs])
        allcalls = allcalls + [call]
        allkernels = allkernels + [comp]

    return allcalls, allkernels

class OPSCCodePrinter(CCodePrinter):
    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return str(float(p)/float(q))

def ccode(expr):
    return OPSCCodePrinter().doprint(expr)
