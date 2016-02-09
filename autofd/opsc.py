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
        lh = flatten(list(list(eq.lhs.atoms(Indexed)) for eq in evals))
        rh = flatten(list(list(eq.rhs.atoms(Indexed)) for eq in evals))
        tot_indexed = list(set(lh+rh))
        libs = set(list(i.base.label for i in lh))
        ribs = set(list(i.base.label for i in rh))
        inouts = libs.intersection(ribs)
        ins = ribs.difference(inouts)
        outs = libs.difference(inouts)
        inouts = list(inouts)
        ins = list(ins)
        outs = list(outs)
        tot_base = ins+outs+inouts
        symbs = flatten(list(list(eq.atoms(Symbol)) for eq in evals))
        Idxs = list(set(list(i.indices for i in (lh+rh))))
        Idxs = flatten(list(i) for i in Idxs)
        labes = list(set(list(i.label for i in Idxs)))
        for i in Idxs:
            labes = labes + [i.lower, i.upper]
        symbs = list(set(symbs).difference(set(labes)).difference(set(inp.constants)).difference(set(tot_base)))
        symbs = list(set(symbs))
        out = []
        symdec = []
        evdict = equations_to_dict(evals)
        for sym in symbs:
            symdec = symdec + ['double %s;' % sym]
            if evdict.get(sym):
                pass
            else:
                raise ValueError("I don't know the formula for %s" % sym)

        # Grid range
        lower = []
        upper = []
        for dim in range(inp.ndim):
            lower = lower + [inp.blockdims[dim].lower - inp.halos[dim]]
            upper = upper + [inp.blockdims[dim].upper + inp.halos[dim]+1]

        for ev in evals:
            code = ccode(ev)
            code = code.replace('==', '=') + END_OF_STATEMENT_DELIMITER['OPSC']
            out = out + [code]
        kercall = []
        kerheader = []
        kernel = []
        kername = inp.kername % inp.kernel_ind
        inp.kernel_ind = inp.kernel_ind + 1
        kerca = 'ops_par_loop(%s, \"%s\", %s[%s], %d, %%s' % (kername, kername, inp.block_name, inp.block, inp.ndim)
        head = 'void %s(' % kername
        if tot_base:
            for ind, v in enumerate(tot_base):
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
                varib = flatten(list(v1 for v1 in tot_indexed if v1.base.label == v))
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

                    sten = ','.join(list(str(t) for t in temp))
                    if inp.stencils.get(sten):
                        sten_na = inp.stencils.get(sten)
                    else:
                        sten_na = inp.sten_name % inp.sten_ind
                        inp.stencils[sten] = sten_na
                        inp.sten_ind = inp.sten_ind + 1

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
                    call = arg_call % ('ops_arg_dat', v, sten_na, datatype, opstype)
                    kercall = kercall + [call]
                    kerheader = kerheader + [headty % v]
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
                    kercall = kercall + [call]
                    kerheader = kerheader + [headty % v]
            iter_range = []
            for dim in range(inp.ndim):
                iter_range = iter_range + [str(lower[dim])] + [str(upper[dim])]
            iter_range = ','.join(iter_range)
            kercall.insert(0, kerca % 'iter_range%d' % inp.iterrange)
            for indno in range(len(kercall)-1):
                kercall[indno] = kercall[indno] + ','
            kercall[-1] = kercall[-1] + ');'
            kercall = ['int iter_range%d[] = {%s};\n' % (inp.iterrange, iter_range)] + kercall
            inp.iterrange = inp.iterrange + 1
            kerheader = head + ', '.join(kerheader) + '){'
            kernel = [kerheader] + symdec + out + ['}']
        else:
            LOG.debug(tot_base)
            pass
        return kercall, kernel
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
    # pprint('\n kernel is')
    # print('\n\n'.join(allkernels))
    # pprint('\n Call is')
    # print('\n\n'.join(allcalls))
    return allcalls, allkernels

class OPSCCodePrinter(CCodePrinter):
    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return str(float(p)/float(q))

def ccode(expr):
    return OPSCCodePrinter().doprint(expr)
