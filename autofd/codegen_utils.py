#!/usr/bin/env python

""" This contains utility routines for code generation. """

from sympy import *
import textwrap

import logging
LOG = logging.getLogger(__name__)

# Language-specific delimiters
COMMENT_DELIMITER = {"OPSC": "//", "F90": "!"}
END_OF_STATEMENT_DELIMITER = {"OPSC": ";", "F90": ""}


def header_code(inp, alg):
    """ Generate header code, containing include/import statements and main program routine, in the language of choice. """

    out = []
    lang = alg.language
    if alg.language == 'OPSC':
        out = out + ['#include <stdlib.h>']
        out = out + ['#include <string.h>']
        out = out + ['#include <math.h>']

        out.append('%s Global Constants in the equations are' % COMMENT_DELIMITER[lang])
        doubles = []
        ints = []
        for con in inp.constants:
            if isinstance(con, Symbol):
                if con.is_integer:
                    ints = ints + ['%s' % con]
                else:
                    doubles = doubles + ['%s' % con]
            elif isinstance(con, Indexed):
                tot = 0
                for ind in con.shape:
                    tot = tot + ind
                if con.is_integer:
                    ints = ints + ['%s[%d]' % (con.base, tot)]
                else:
                    doubles = doubles + ['%s[%d]' % (con.base, tot)]
        if ints:
            out = out + ['int %s %s' % (', '.join(ints), END_OF_STATEMENT_DELIMITER[lang])]
        if doubles:
            out = out + ['double %s %s' % (', '.join(doubles), END_OF_STATEMENT_DELIMITER[lang])]

        out = out + ['// OPS header file']
        out = out + ['#define OPS_%sD' % inp.ndim]
        out = out + ['#include "ops_seq.h"']
        out = out + ['#include "auto_kernel.h"']
        out = out + ['%s main program start' % COMMENT_DELIMITER[lang]]
        out = out + ['int main (int argc, char **argv) {']
    elif lang == 'F90':
        doubles = []
        ints = []
        for con in inp.constants:
            if isinstance(con, Symbol):
                if con.is_integer:
                    ints = ints + ['%s' % con]
                else:
                    doubles = doubles + ['%s' % con]
            elif isinstance(con, Indexed):
                tot = 0
                LOG.debug(fcode(con))
                for ind in con.shape:
                    tot = tot + ind
                if con.is_integer:
                    ints = ints + ['%s(%d)' % (con.base, tot)]
                else:
                    doubles = doubles + ['%s(%d)' % (con.base, tot)]
        if ints:
            out = out + ['integer :: %s %s' % (', '.join(ints), END_OF_STATEMENT_DELIMITER[lang])]
        if doubles:
            out = out + ['real(8) :: %s %s' % (', '.join(doubles), END_OF_STATEMENT_DELIMITER[lang])]
        inp.module.append('\n'.join(out))
        out = []
        # TODO spj change module name later as of now using the same module name parammod
        out = out + ['program codegen']
        out = out + ['use param_mod']
        out = out + ['IMPLICIT NONE']

    else:
        raise ValueError('Implement header declaration for the %s language' % lang)

    return out


def loop(indices, alg):
    """ Generate a 'for' loop in the desired output language.
    
    :arg indices: the loop indices.
    :returns: headers and footers of the loops.
    """

    header = []
    footer = []

    if alg.language == 'OPSC':
        comment = '%s loop start %s' % (COMMENT_DELIMITER[alg.language], indices[0])
        header.append(comment)
        for dim in range(0, len(indices)):
            temp = indices[dim]
            header.append('for(int %s=%s; %s<%s; %s++){' % (temp, temp.lower, temp, temp.upper, temp))
            footer.append('}')
        comment = '%s loop end %s' % (COMMENT_DELIMITER[alg.language], indices[0])
        footer.append(comment)

    elif alg.language == 'F90':
        comment = '%s loop start %s' % (COMMENT_DELIMITER[alg.language], indices[0])
        header.append(comment)
        for dim in reversed(range(0, len(indices))):
            temp = indices[dim]
            header.append('do %s=%s,%s' % (temp, temp.lower, temp.upper))
            footer.append('enddo')
        comment = '%s loop end %s' % (COMMENT_DELIMITER[alg.language], indices[0])
        footer.append(comment)

    else:
        raise ValueError('Implement loop declaration for the %s language.' % alg.language)

    return header, footer


def defdec(inp, alg):
    """ Define and declare variables. """

    out = []
    lang = alg.language

    if alg.language == 'OPSC':
        # Define inputs to the code
        block_name = inp.block_name
        totblock = inp.block.upper - inp.block.lower
        ind = inp.blockdims
        inputs = []
        inputs = inputs + ['int %s = %d %s' % (totblock, inp.nblocks, END_OF_STATEMENT_DELIMITER[lang])]

        # change this
        inputs = inputs + ['int %s %s' % (','.join(list('%s[%s]' % ('nx%dp' % dim, totblock) for dim in range(inp.ndim))), END_OF_STATEMENT_DELIMITER[lang])]

        if not inp.multiblock:
            inputs = inputs + ['int %s = %d%s' % (inp.block, inp.nblocks-1, END_OF_STATEMENT_DELIMITER[lang])]
            inputs = inputs + ['\n'.join(list('%s = ;' % inp.grid[dim+1] for dim in range(inp.ndim)))]
        else:
            inputs = inputs + ['%s Write the block dimensions here' % (COMMENT_DELIMITER[lang])]
            inputs = inputs + ['\n\n']
            inputs = inputs + ['%s Writing the block dimensions ends here' % (COMMENT_DELIMITER[lang])]
        # inputs

        # Declare Constants in OPS format
        out = out + ['ops_init(argc,argv,1)%s' % END_OF_STATEMENT_DELIMITER[lang]]
        for constant in inp.constants:
            if isinstance(constant, Symbol):
                inputs = inputs + ['%s = %s' % (constant, END_OF_STATEMENT_DELIMITER[lang])]
            elif isinstance(constant, Indexed):
                tot = 0
                for inde in constant.shape:
                    tot = tot + inde
                for no in range(tot-1):
                    inputs = inputs + ['%s[%d] = %s' % (constant.base, no, END_OF_STATEMENT_DELIMITER[lang])]
        for con in inp.constants:
            if con.is_Symbol:
                if con.is_integer:
                    dtype = 'int'
                else:
                    dtype = 'double'
                out = out + ['ops_decl_const(\"%s\" , 1, \"%s\", &%s)%s' % (con, dtype, con, END_OF_STATEMENT_DELIMITER[lang])]

        # Declare block
        out.append('ops_block *%s = (ops_block *)malloc(%s*sizeof(ops_block*));' % (block_name, totblock))
        out.append('char buf[100];')
        if inp.multiblock:
            stloop, enloop = loop([inp.block], alg)
        else:
            stloop = ['\n']
            enloop = ['\n']
        out = out + stloop
        out.append('sprintf(buf,\"%s[%%d]\",%s);' % (block_name, inp.block))
        out.append('%s[%s] = ops_decl_block(%d,buf);' % (block_name, inp.block, inp.ndim))

        out = out + enloop

        # Declare data files
        for da in inp.dats:
            out.append('ops_dat *%s = (ops_dat *)malloc(%s*sizeof(ops_dat*));' % (da.base, totblock))
        out = out + ['int d_p[%s]   = {%s};' % (inp.ndim, ','.join(list('%d' % (inp.halos[dim]) for dim in range(inp.ndim))))]
        out = out + ['int d_m[%d]   = {%s}; ' % (inp.ndim, ','.join(list('-%d' % (inp.halos[dim]) for dim in range(inp.ndim))))]
        out = out + ['int base[%d]  = {%s};' % (inp.ndim, ','.join(list('%d' % (0) for dim in range(inp.ndim))))]
        out = out + ['double* temp = NULL;']
        out = out + stloop

        size = ','.join(list('%s' % (ind[dim].upper - ind[dim].lower + 1) for dim in range(inp.ndim)))
        out = out + ['int size[%d] = {%s};' % (inp.ndim, size)]
        for da in inp.dats:
            out.append('sprintf(buf,\"%s[%%d]\",%s);' % (da.base, inp.block))
            out = out + ['%s[%s] = ops_decl_dat(%s[%s], 1, size, base, d_m, d_p, temp, "double", buf);' % (da.base, inp.block, block_name, inp.block)]
        out = out + enloop

        for key, value in inp.stencils.iteritems():
            npts = divmod(len(key.split(',')), inp.ndim)
            if npts[1] == 0:
                pass
            else:
                raise ValueError('number of points are not a multiple of dimensions')
            tname = 'sten_%s' % value
            out = out + ['int %s[] = {%s}%s' % (tname, key, END_OF_STATEMENT_DELIMITER[lang])]
            sname = 'ops_stencil %s = ops_decl_stencil(%d, %d, %s, \"%s\")%s'
            out = out + [sname % (value, inp.ndim, npts[0], tname, tname, END_OF_STATEMENT_DELIMITER[lang])]
        out = out + ['\n\n']+inp.bcdecl + ['\n\n'] + ['ops_partition("");']
        out = inputs + out

    elif lang == 'F90':
        inputs = []
        # These are the grid dimensions etc..
        inputs = inputs + ['integer :: %s %s' % (','.join(list('%s' % ('nx%dp' % dim) for dim in range(inp.ndim))), END_OF_STATEMENT_DELIMITER[lang])]
        inputs = inputs + ['\n'.join(list('%s = ' % inp.grid[dim+1] for dim in range(inp.ndim)))]
        for con in inp.constants:
            if isinstance(con, Symbol):
                inputs = inputs + ['%s = %s' % (con, END_OF_STATEMENT_DELIMITER[lang])]
            elif isinstance(con, Indexed):
                tot = 0
                for inde in con.shape:
                    tot = tot + inde
                for no in range(tot):
                    inputs = inputs + ['%s(%d) = %s' % (con.base, no, END_OF_STATEMENT_DELIMITER[lang])]
        out = []
        dimen = ','.join(list(':' for dim in range(inp.ndim)))

        for da in inp.dats:
            out.append('%s' % (da.base))
        out = 'real(8), allocatable, dimension (%s) :: ' % (dimen) + ', '.join(out)
        out = ' &\n'.join(textwrap.wrap(out, width=70, break_long_words=False))
        inp.module.append(out)
        out = []
        # Allocate stuff
        ind = inp.blockdims
        sz = []
        for dim in range(inp.ndim):
            l = ind[dim].lower - inp.halos[dim]
            u = ind[dim].upper + inp.halos[dim]
            s = str(l) + ':' + str(u)
            sz.append(s)
        sz = ', '.join(sz)
        for da in inp.dats:
            out = out + ['allocate (%s(%s))' % (da.base, sz)]
        out = inputs + out
        out = list(' '*6+ou for ou in out)
        # out = '\n'.join(out)
        # print(out)
        # TODO write this to take care of long allocations

    else:
        raise ValueError('Implement %s in the Definitions and declaration' % alg.language)

    return out


def footer_code(inp, alg):
    out = []
    if alg.language == 'OPSC':
        out = out + ['ops_printf(\" finished running the code\\n\");', 'ops_exit();', '}']
    elif alg.language == 'F90':
        out = out + ['end program']
    else:
        raise ValueError('Implement %s in the footer code' % alg.language)
    return out


def indent_code(code):
    """ Indent the code.
    
    :arg code: a string of code or a list of code lines.
    :returns: the indented code
    """

    if isinstance(code, string_types):
        code_lines = indent_code(code.splitlines(True))
        return ''.join(code_lines)

    tab = "   "
    inc_token = ('{', '(', '{\n', '(\n')
    dec_token = ('}', ')')

    code = [line.lstrip(' \t') for line in code]

    increase = [int(any(map(line.endswith, inc_token))) for line in code]
    decrease = [int(any(map(line.startswith, dec_token)))
                for line in code]

    pretty = []
    level = 0
    for n, line in enumerate(code):
        if line == '' or line == '\n':
            pretty.append(line)
            continue
            level -= decrease[n]
            pretty.append("%s%s" % (tab*level, line))
            level += increase[n]
    return pretty


def bc_periodic_OPSC():

    return


def bcs(inp, alg):
    """ Generate code for the boundary conditions. """

    inp.bcdecl = []
    inp.bccall = []
    inp.bcexind = 0
    lang = alg.language
    if len(alg.bcs) == inp.ndim:
        pass
    elif (len(alg.bcs) > inp.ndim):
        raise ValueError('There are more boundary conditions than the number of dimensions')
    elif (len(alg.bcs) < inp.ndim):
        raise ValueError('There are less boundary conditions than the number of dimensions')
    inp.bc_appl = []
    for dim in range(inp.ndim):
        bc = alg.bcs[dim]
        if bc[0] == bc[1] and bc[0] == 'periodic' and alg.language == 'OPSC' and not inp.multiblock:
            out = []
            iter_size = list(te for te in inp.gridhalo)

            LOG.debug('periodic bc in x%d- direction' % dim)
            halo = inp.halos[dim]  # get the number of halo points
            iter_size[dim] = halo  # the no of halos to be transfered
            from_base = list(-halo for te in range(inp.ndim))
            to_base = list(-halo for te in range(inp.ndim))

            # Now copy from data at the end of domain to the first
            l = inp.blockdims[dim].lower - halo
            u = inp.blockdims[dim].upper - halo + 1
            from_base[dim] = u
            to_base[dim] = l
            iter_size = ','.join(list(str(te) for te in iter_size))
            fro = ','.join(list(str(te) for te in from_base))
            to = ','.join(list(str(te) for te in to_base))
            # dire = ','.join(['1','2'])
            dire = ','.join(str(i) for i in range(1, inp.ndim+1))
            halo_count = 0
            stloop, enloop = loop([inp.block], alg)
            out = out + stloop
            out = out + ['int off = 0;']
            out = out + ['int halo_iter[] = {%s}%s' % (iter_size, END_OF_STATEMENT_DELIMITER[alg.language])]
            out = out + ['int from_base[] = {%s}%s' % (fro, END_OF_STATEMENT_DELIMITER[alg.language])]
            out = out + ['int to_base[] = {%s}%s' % (to, END_OF_STATEMENT_DELIMITER[alg.language])]
            out = out + ['int dir[] = {%s}%s' % (dire, END_OF_STATEMENT_DELIMITER[alg.language])]
            haloname = 'halos_x%d_%s' % (dim, bc[0])
            # halos[off++] = ops_decl_halo(u[i-1+ngrid_x*j], u[i+ngrid_x*j], halo_iter, base_from, base_to, dir, dir);
            halo_decl = '%s[off++] = ops_decl_halo(%%s[%%s],%%s[%%s],halo_iter, from_base, to_base, dir, dir);' % haloname
            for con in inp.conser:
                out = out + [halo_decl % (con.base.label, inp.block, con.base.label, inp.block)]
                halo_count = halo_count + 1

            # Copy the first n terms to the Rightend halos
            l = inp.blockdims[dim].lower
            u = inp.blockdims[dim].upper + 1
            from_base[dim] = l
            to_base[dim] = u
            fro = ','.join(list(str(te) for te in from_base))
            to = ','.join(list(str(te) for te in to_base))
            out = out + ['from_base[%d] = %s%s' % (dim, l, END_OF_STATEMENT_DELIMITER[alg.language])]
            out = out + ['to_base[%d] = %s%s' % (dim, u, END_OF_STATEMENT_DELIMITER[alg.language])]
            for con in inp.conser:
                out = out + [halo_decl % (con.base.label, inp.block, con.base.label, inp.block)]
                halo_count = halo_count + 1
            out = ['ops_halo *%s = (ops_halo *)malloc(%d*sizeof(ops_halo *));' % (haloname, halo_count)] + out
            out = out + enloop
            out = out + ['ops_halo_group grp_%s = ops_decl_halo_group(%s,%s);' % (haloname, halo_count, haloname)]
            inp.bcdecl.append([out])
            inp.bccall.append('ops_halo_transfer(grp_%s);' % haloname)

        else:
            for loc_bound, bound in enumerate(bc):
                out = []
                halo_count = 0
                boundary = 'Bc_x%d_%%s_%s' % (dim, bound)
                iter_size = list(te for te in inp.gridhalo)
                halo = inp.halos[dim]
                from_base = list(-halo for te in range(inp.ndim))
                to_base = list(-halo for te in range(inp.ndim))
                if loc_bound == 0:
                    bcloc = boundary % 'min'
                    val = inp.blockdims[dim].lower
                    dire = -1
                elif loc_bound == 1:
                    bcloc = boundary % 'max'
                    val = inp.blockdims[dim].upper
                    dire = 1
                else:
                    raise ValueError('undefined bc')
                stloop, enloop = loop([inp.block], alg)
                out = out + ['%s\n%s Boundary condition %s\n%s' % (COMMENT_DELIMITER[lang], COMMENT_DELIMITER[lang], bcloc, COMMENT_DELIMITER[lang])]
                out = out + stloop
                te1 = ','.join(list(str(te) for te in iter_size))
                out = out + ['int off = 0;']
                out = out + ['int halo_iter[] = {%s}%s' % (te1, END_OF_STATEMENT_DELIMITER[lang])]
                te1 = ','.join(list(str(te) for te in from_base))
                out = out + ['int from_base[] = {%s}%s' % (te1, END_OF_STATEMENT_DELIMITER[lang])]
                te1 = ','.join(list(str(te) for te in to_base))
                out = out + ['int to_base[] = {%s}%s' % (te1, END_OF_STATEMENT_DELIMITER[lang])]
                out = out + ['int dir[] = {%s}%s' % (','.join(str(i) for i in range(1, inp.ndim+1)), END_OF_STATEMENT_DELIMITER[lang])]
                halo_decl = '%s[off++] = ops_decl_halo(%%s[%%s],%%s[%%s],halo_iter, from_base, to_base, dir, dir);' % bcloc
                LOG.debug(bcloc)

                if bound == 'symmetry':
                    raise ValueError('Implementation %s bc' % bound)
                elif bound == 'zero_grad':
                    # Change the values depending on the direction
                    iter_size[dim] = 1
                    structu = '%s[%d] = %s%s'
                    out = out + [structu % ('halo_iter', dim, iter_size[dim], END_OF_STATEMENT_DELIMITER[lang])]
                    LOG.debug(bound, halo, val)
                    LOG.debug('from is', val - dire*1)
                    for d in range(0, halo+1):
                        from_base[dim] = val - dire*1
                        to_base[dim] = val + dire*d
                        out = out + [structu % ('from_base', dim, from_base[dim], END_OF_STATEMENT_DELIMITER[lang])]
                        out = out + [structu % ('to_base', dim, to_base[dim], END_OF_STATEMENT_DELIMITER[lang])]
                        for con in inp.conser:
                            out = out + [halo_decl % (con.base.label, inp.block, con.base.label, inp.block)]
                            halo_count = halo_count + 1
                else:
                    raise ValueError("Don't know implementation of the boundary condition of type %s." % bound)

                out = ['ops_halo *%s = (ops_halo *)malloc(%d*sizeof(ops_halo *));' % (bcloc, halo_count)] + out
                out = out + enloop
                out = out + ['ops_halo_group grp_%s = ops_decl_halo_group(%s,%s);' % (bcloc, halo_count, bcloc)]
                inp.bcdecl.append([out])
                inp.bccall.append('ops_halo_transfer(grp_%s);' % bcloc)

    return


def after_time(inp, alg):
    out = []
    out = out + write_state(inp, alg)
    return out

def write_state(inp, alg):
    """ Write the current state of the simulation (i.e. all of the fields/conservative variables) to disk in HDF5 format. """
    out = []
    if alg.language == 'OPSC':
        # First write out the block(s)
        block_to_hdf5 = "ops_fetch_block_hdf5_file(%s[%s], \"state.h5\");" % (inp.block_name, inp.block)
        out.append(block_to_hdf5)

        # Then write out each field.
        for c in inp.conser:
            convervative_variable_to_hdf5 = "ops_fetch_dat_hdf5_file(%s[%s], \"state.h5\");" % (c.base.label, inp.block)
            out.append(convervative_variable_to_hdf5)

    elif alg.language == 'F90':
        #TODO: Write HDF5 output for Fortran 90.
        out = out + ['%s Write data file output here' % COMMENT_DELIMITER[alg.language]]
    else:
        raise ValueError('Implement output writing for language %s' % alg.language)

    return out


def write_final_code(template, codes, main, routines, lang):
    import collections
    od = collections.OrderedDict(sorted(template.items()))
    main_code = []
    for key, val in od.iteritems():
        # Go ahead with writing code
        if codes.get(val):
            va = codes.get(val)
            va = flatten(va)
            main_code = main_code + va + ['\n']
        else:
            main_code = main_code + ['%s There is no code provided for %s part in the algorithm\n' % (COMMENT_DELIMITER[lang], val)] + ['\n']

    main_code = flatten(main_code)
    main.write('\n'.join(main_code))

    if lang == 'OPSC':
        temp = ['#ifndef kernels_KERNEL_H \n#define kernels_KERNEL_H\n']  # + codes['kernels'] + ['\n#endif']
        routines.write('\n'.join(temp))
        temp = flatten(codes['kernels'])
        temp = flatten(temp)
        routines.write('\n'.join(temp))
        routines.write('\n#endif')
    elif lang == 'F90':
        temp = flatten(codes['kernels'])
        routines.write('\n'.join(temp))
    else:
        raise ValueError('Implement %s in write_final_code' % lang)

    return


def write_module(mods, modfile):
    temp = []
    pprint(mods)
    temp = temp + ['module param_mod'] + mods + ['end module param_mod']
    modfile.write('\n'.join(temp))

    return

