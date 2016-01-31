#!/usr/bin/env python

from sympy import *
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_application)
transformations = standard_transformations + (implicit_application,)

# AutoFD functions
from .codegen_utils import loop, COMMENT_DELIMITER, END_OF_STATEMENT_DELIMITER

import logging
LOG = logging.getLogger(__name__)


def F90_write_kernel(eqs, inp, alg):
    def get_kernel(evals):
        lh = flatten(list(list(eq.lhs.atoms(Indexed)) for eq in evals))
        rh = flatten(list(list(eq.rhs.atoms(Indexed)) for eq in evals))
        tot_indexed = list(set(lh+rh))
        libs = set(list(i.base.label for i in lh))
        ribs = set(list(i.base.label for i in rh))
        indi = list(set(list(i.indices for i in (lh+rh))))
        inouts = libs.intersection(ribs)
        ins = ribs.difference(inouts)
        outs = libs.difference(inouts)
        inouts = list(inouts)
        ins = list(ins)
        outs = list(outs)
        tot_base = ins+outs+inouts
        if len(indi) == 1:
            pass
        out = []
        # grid range
        lower = []
        upper = []
        for dim in range(inp.ndim):
            lower = lower + [inp.blockdims[dim].lower - inp.halos[dim]]
            upper = upper + [inp.blockdims[dim].upper + inp.halos[dim]]
        # LOG.debug(lower,upper)
        for ev in evals:
            code = fcode(ev, source_format='free', standard=95)
            code = code.replace('==', '=') + END_OF_STATEMENT_DELIMITER['F90']
            out = out + [code]
        kercall = []
        kerheader = []
        kernel = []
        kername = inp.kername % inp.kernel_ind
        inp.kernel_ind = inp.kernel_ind + 1
        inpargs = []
        # kerheader = kerheader+ [''*78]+ [''*14subroutine %s'%kername]
        if tot_base:
            for ind, v in enumerate(tot_base):
                varib = flatten(list(v1 for v1 in tot_indexed if v1.base.label == v))
                varib = list(set(varib))
                variabind = flatten(list(v1.indices) for v1 in varib)
                variabind = list(set(variabind))

                if all(va.upper == inp.block.upper for va in variabind):
                    indexes = list(va for va in variabind)
                    for dim in range(inp.ndim):
                        indexes = list(str(te).replace('x%d' % dim, 'i%d' % dim) for te in indexes)
                    for inde in range(len(indexes)):
                        for ou in range(len(out)):
                            if isinstance(indexes[inde], tuple):
                                new = '%s' % (indexes[inde])
                            else:
                                new = '%s' % (indexes[inde])
                            old = ('%s' % (variabind[inde]))
                            # out[ou] = out[ou].replace(old, new)
                    for dim in range(inp.ndim):
                        indexes = list(str(te).replace('i%d' % dim, '0') for te in indexes)
                    indexes = list(parse_expr(v) for v in indexes)
                    indexes = indexes + [parse_expr(', '.join(list(str(0) for dim in range(inp.ndim))))]
                    indexes = list(set(indexes))
                    if inp.ndim > 1:
                        for dim in range(inp.ndim):
                            indexes = sorted(indexes, key=lambda indexes: indexes[dim])
                    else:
                        indexes = [sorted(indexes)]
                    # update range on which the loop to be iterated
                    if len(indexes) == 1:
                        for dim in range(inp.ndim):
                            lower[dim] = lower[dim] - indexes[0][dim]
                            upper[dim] = upper[dim] - indexes[0][dim]
                    else:
                        for dim in range(inp.ndim):
                            lower[dim] = lower[dim] - indexes[0][dim]
                            upper[dim] = upper[dim] - indexes[-1][dim]
                else:
                    indexes = list(va for va in variabind)
                    inpargs = inpargs + indexes
            inpargs = list(set(inpargs))
            n_inden = 6
            if inpargs:
                temp = ','.join(list(str(t.label) for t in inpargs))
                kerheader = [kerheader] + ['subroutine %s(%s)' % (kername, temp)] + ['use param_mod'] + ['IMPLICIT NONE'] + ['\n']
                kercall = kercall + [' '*n_inden + 'call %s(%s)' % (kername, temp)]
                ints = []
                doubs = []
                for t in inpargs:
                    if t.is_integer:
                        ints = ints + [str(t.label)]
                    else:
                        doubs = doubs + [str(t.label)]
                if ints:
                    kerheader = [kerheader] + ['integer :: %s' % (','.join(ints))]
                if doubs:
                    kerheader = [kerheader] + ['real(8) :: %s' % (','.join(doubs))]
            else:
                kerheader = [kerheader] + ['subroutine %s' % (kername)] + ['use param_mod'] + ['IMPLICIT NONE'] + ['\n']
                kercall = kercall + [' '*n_inden + 'call %s' % (kername)]

            # if len(indi) - len(inpargs) == 1:
                # inde = ','.join(list('i%d'%dim for dim in range(inp.ndim)))
                # new = ','.join(list(':' for dim in range(inp.ndim)))
                # for ou in range(len(out)):
                    # out[ou] = out[ou].replace(inde,new)
                # kerheader =[kerheader]
            # else:
            kerheader = [kerheader] + ['integer :: %s' % (','.join(list('x%d' % dim for dim in range(inp.ndim))))]
            # kerheader =[kerheader] +  ['IMPLICIT NONE'] +['\n']
            newdims = []
            for dim in range(inp.ndim):
                newdims = newdims + [Idx('x%d' % dim, (lower[dim], upper[dim]))]
            st, en = loop(newdims, alg)
            out = st + out + en

            kernel = ['%s Start  sub routine \n' % COMMENT_DELIMITER['F90']]
            kernel = kernel + kerheader + out + ['end subroutine']
            # kernel = '\n'.join(kernel)
            kernel = kernel + ['%s End sub routine \n\n\n' % COMMENT_DELIMITER['F90']]
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
    return allcalls, allkernels
