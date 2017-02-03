#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

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

from sympy.printing.ccode import CCodePrinter
import os
import logging
LOG = logging.getLogger(__name__)
BUILD_DIR = os.getcwd()


class OPSCCodePrinter(CCodePrinter):

    """ Prints OPSC code. """

    def __init__(self, Indexed_accs, constants):
        """ Initialise the code printer. """

        settings = {}
        CCodePrinter.__init__(self, settings)

        # Indexed access numbers are required in dictionary
        self.Indexed_accs = Indexed_accs
        self.constants = constants

    def _print_ReductionVariable(self, expr):
        return '*%s' % str(expr)

    def _print_Rational(self, expr):
        if self.constants is not None:
            if expr in self.constants.keys():
                return self.constants[expr]
            else:
                con = 'rc%d' % len(self.constants.keys())
                self.constants[expr] = con
                return self.constants[expr]
        else:
            p, q = int(expr.p), int(expr.q)
            return '%d.0/%d.0' % (p, q)

    def _print_Mod(self, expr):
        args = map(ccode, expr.args)
        args = [x for x in args]
        result = ','.join(args)
        result = 'fmod(%s)' % result
        return result

    def _print_Max(self,expr):
        nargs = len(expr.args)
        args_code = [self._print(a) for a in expr.args]
        if nargs == 2:
            args_code = ', '.join(args_code)
            return "MAX(%s)"%(args_code)
        # Need to come up with a better IDEA For this FIXME
        elif nargs == 3:
            a = args_code[0]; c = args_code[2];
            b = args_code[1];
            return "MAX(MAX(%s,%s), MAX(%s,%s))"%(a,b,b,c)
        elif nargs == 4:
            a = args_code[0]; c = args_code[2];
            b = args_code[1]; d = args_code[3];
            return "MAX(MAX(MAX(%s,%s), MAX(%s,%s)), MAX(%s,%s))"%(a,b,b,c,c,d)
        elif nargs == 5:
            a = args_code[0]; c = args_code[2]; e = args_code[4]
            b = args_code[1]; d = args_code[3];
            return "MAX(MAX(MAX(%s,%s), MAX(%s,%s)), MAX(MAX(%s,%s),MAX(%s,%s)))"%(a,b,b,c,c,d,d,e)
        elif nargs == 6:
            a = args_code[0]; c = args_code[2]; e = args_code[4]
            b = args_code[1]; d = args_code[3]; f = args_code[5]
            return "MAX(MAX(MAX(%s,%s), MAX(%s,%s)), MAX(MAX(%s,%s),MAX(%s,%s)))"%(a,b,b,c,c,d,d,e)
        else:
            raise ValueError("Max for arguments %d is not defined in code printer or OPS"%nargs)

    def _print_Indexed(self, expr):
        """ Print out an Indexed object.

        :arg expr: The Indexed expression.
        :returns: The indexed expression, as OPSC code.
        :rtype: str
        """

        # Replace the symbols in the indices that are not time with `zero'
        indices = [ind for ind in expr.indices if ind != EinsteinTerm('t')]
        for number, index in enumerate(indices):
            for sym in index.atoms(Symbol):
                indices[number] = indices[number].subs({sym: 0})
        if self.Indexed_accs:
            if self.Indexed_accs[expr.base]:
                out = "%s[%s(%s)]" % (self._print(expr.base.label), self.Indexed_accs[expr.base], ','.join([self._print(index) for index in indices]))
            else:
                out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(index) for index in indices]))
        else:
            out = "%s[%s]" % (self._print(expr.base.label), ','.join([self._print(index) for index in indices]))
        return out


def pow_to_constant(expr, constants):
    from sympy.core.function import _coeff_isneg
    # Only negative powers i.e they correspond to division and they are stored into constant arrays
    inverse_terms = {}
    for at in expr.atoms(Pow):
        if _coeff_isneg(at.exp) and not at.base.atoms(Indexed) and not at.base.atoms(GridVariable):
            if not at in constants.keys():
                constants[at] = 'rinv%d' % len(constants.keys())
            inverse_terms[at] = constants[at]
    expr = expr.subs(inverse_terms)
    return expr, constants


def ccode(expr, Indexed_accs=None, constants=None):
    """ Create an OPSC code printer object and write out the expression as an OPSC code string.

    :arg expr: The expression to translate into OPSC code.
    :arg Indexed_accs: Indexed OPS_ACC accesses.
    :arg constants: Constants that should be defined at the top of the OPSC code.
    :returns: The expression in OPSC code.
    :rtype: str
    """
    if isinstance(expr, Eq):
        if constants:
            expr, constants = pow_to_constant(expr, constants)
        code_print = OPSCCodePrinter(Indexed_accs, constants)
        code = code_print.doprint(expr.lhs) \
            + ' = ' + OPSCCodePrinter(Indexed_accs, constants).doprint(expr.rhs)
        return code, code_print.constants
    return OPSCCodePrinter(Indexed_accs, constants).doprint(expr)

from opensbli.core.kernel import Kernel
from opensbli.core.algorithm import Loop
from opensbli.core.opensbliobjects import *

class BeforeTimeOpsc():
    def __init__(self):
        self.components = []
        return
    def add_components(self):
        if isinstance(components, list):
            self.components += components
        else:
            self.components += [components]
class DeclareDataset(object):
    def __init__(self, dataset, blocknumber):
        self.dataset = dataset
        self.block_number = blocknumber
        return
class WriteString(object):
    def __init__(self, string):
        if isinstance(string, list):
            self.components = string
        elif isinstance(string, str):
            self.components = [string]
        else:
            raise ValueError("")
        return
    def __str__(self):
        return '\n'.join(self.components)
    def _print(self):
        return str(self)
    @property
    def opsc_code(self):
        return ['\n'.join(self.components)]


class OPSC(object):

    def __init__(self, algorithm):
        """ Generating an OPSC code from the algorithm"""
        if not algorithm.MultiBlock:
            self.MultiBlock = False
            self.datasets_to_declare = []
            self.Rational_constants = set()
            self.constants = set()
            self.stencils_to_declare = set()
            self.dtype = algorithm.dtype
            #self.generate_OPSC_dependants_sb(algorithm)
            def_decs = self.opsc_def_decs(algorithm)
            algorithm.prg.components = def_decs + algorithm.prg.components
            code = algorithm.prg.opsc_code
            f = open('test.cpp', 'w')
            f.write('\n'.join(code))
            f.close()
        return

    def opsc_def_decs(self, algorithm):
        defs = []
        decls = []
        # Add OPS_init to the declarations as it should be called before all ops
        decls += self.ops_init()
        # First process all the constants in the definitions
        for d in algorithm.defnitionsdeclarations.components:
            if isinstance(d, Constant):
                defs += [self.define_constants(d)]
                decls += [self.declare_ops_constants(d)]
        # Once the constants are done define and declare OPS dats
        output = defs + decls
        defs = []
        decls = []
        # Define and declare blocks
        for b in algorithm.block_descriptions:
            output += self.define_block(b)
            output += self.declare_block(b)
        #output += defs + decls
        # Define and declare datasets on each block
        for d in algorithm.defnitionsdeclarations.components:
            if isinstance(d, DataSetBase):
                #output += self.define_dataset(d)
                output += self.declare_dataset(d)
        #output = output + defs + decls
        # Define stencils
        output += [WriteString("// Define and declare stencils")]
        # Loop through algorithm components to include any halo exchanges
        from opensbli.core.bcs import Exchange
        exchange_list = self.loop_alg(algorithm, Exchange)
        for e in exchange_list:
            call, code = self.bc_exchange_call_code(e)
            output += [code]
        return output

    def Exchange_code(self, e):
        #out =
        return
    def bc_exchange_call_code(self, instance):
        off = 0
        halo = 'halo'
        instance.transfer_size = instance.transfer_from
        # Name of the halo exchange
        name = instance.name
        #self.halo_exchange_number = self.halo_exchange_number + 1
        code = ['// Boundary condition exchange code' ]
        code += ['ops_halo_group %s %s' % (name, ";")]
        code += ["{"]
        code += ['int halo_iter[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_size]), ";")]
        code += ['int from_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_from]), ";")]
        code += ['int to_base[] = {%s}%s' % (', '.join([str(s) for s in instance.transfer_to]), ";")]
        # dir in OPSC. FIXME: Not sure what it is, but 1 to ndim works.
        code += ['int dir[] = {%s}%s' % (', '.join([str(ind+1) for ind in range(len(instance.transfer_to))]), ";")]
        # Process the arrays
        for arr in instance.transfer_arrays:
            code += ['ops_halo %s%d = ops_decl_halo(%s, %s, halo_iter, from_base, to_base, dir, dir)%s'
                     % (halo, off, arr.base, arr.base, ";")]
            off = off+1
        code += ['ops_halo grp[] = {%s}%s' % (','.join([str('%s%s' % (halo, of)) for of in range(off)]), ";")]
        code += ['%s = ops_decl_halo_group(%d,grp)%s' % (name, off, ";")]
        code += ["}"]
        # Finished OPS halo exchange, now get the call
        instance.call_name = 'ops_halo_transfer(%s)%s' % (name, ";")
        call = ['// Boundary condition exchange calls' , 'ops_halo_transfer(%s)%s' % (name, ";")]
        for no, c in enumerate(code):
            code[no] = WriteString(c)
        return call, code

    def loop_alg(self, algorithm, type_of_component):
        type_list = []
        def _generate(components, type_list):
            for component1 in components:
                if isinstance(component1, type_of_component):
                    type_list += [component1]
                elif isinstance(component1, Loop):
                    return _generate(component1.components, type_list)

        _generate(algorithm.prg.components, type_list)
        return type_list

    def ops_init(self, diagnostics_level=None):
        """ The default diagnostics level is 1, which offers no diagnostic information and should be used for production runs.
        Refer to OPS user manual for more information.

        :arg int diagnostics_level: The diagnostics level. If None, the diagnostic level defaults to 1.
        :returns: The call to ops_init.
        :rtype: list
        """
        out = [WriteString('// Initializing OPS ')]
        if diagnostics_level:
            self.ops_diagnostics = True
            return out + [WriteString('ops_init(argc,argv,%d);' % (diagnostics_level))]
        else:
            self.ops_diagnostics = False
            return out + [WriteString('ops_init(argc,argv,1);')]

    def define_block(self, b):
        if not self.MultiBlock:
            return [WriteString("ops_block %s;"%b.block_name)]
        else:
            raise NotImplementedError("")

    def declare_block(self, b):
        if not self.MultiBlock:
            string = '%s = ops_decl_block(%d, \"%s\");' % (b.block_name, b.ndim, b.block_name)
            return [WriteString(string)]
        else:
            raise NotImplementedError("")

    def define_constants(self, c):
        if isinstance(c, ConstantObject):
            if c.value:
                return WriteString("%s=%s;"%(str(c), ccode(c.value)))
            else:
                return WriteString("%s=%s;"%(str(c), "Input"))
        elif isinstance(c, ConstantIndexed):
            if c.value:
                raise NotImplementedError("")
            else:
                raise NotImplementedError("")
        else:
            print c
            raise ValueError("")
    def declare_ops_constants(self, c):
        if isinstance(c, ConstantObject):
            return WriteString("ops_decl_const(\"%s\" , 1, \"%s\", &%s);"%(str(c), c.dtype, str(c)))
        return
    def declare_inline_array(self, dtype, name, values):
        return WriteString('%s %s[] = {%s};' % (dtype, name, ', '.join([str(s) for s in values])))
    def update_inline_array(self, name, values):
        out = []
        for no,v in enumerate(values):
            out += [WriteString("%s[%d] = %s;"%(name, no, value))]
        return out
    def define_dataset(self, dset):
        if not self.MultiBlock:
            return [WriteString("ops_dat %s;" %(dset))]
    def get_max_halos(self, halos):
        halo_m = []
        halo_p = []
        for direction in range(len(halos)):
            max_halo_direction = []
            if halos[direction][0]:
                hal = [d.get_halos(0) for d in halos[direction][0]]
                halo_m += [min(hal)]
            else:
                halo_m += [0]
            if halos[direction][1]:
                hal = [d.get_halos(1) for d in halos[direction][1]]
                halo_p += [max(hal)]
            else:
                halo_p += [0]
        return halo_m, halo_p
    def declare_dataset(self, dset):
        hm, hp = self.get_max_halos(dset.halo_ranges)
        halo_p = self.declare_inline_array("int", "halo_p", hp)
        halo_m = self.declare_inline_array("int", "halo_m", hm)
        sizes = self.declare_inline_array("int", "size", [str(s) for s in (dset.size)])
        base = self.declare_inline_array("int", "base", [0 for i in range(len(dset.size))])
        if dset.dtype:
            dtype = dset.dtype
        else:
            dtype = self.dtype
        value = WriteString("%s value = NULL;"%dtype)
        temp = 'ops_dat %s = ops_decl_dat(%s, 1, size, base, halo_m, halo_p, value, \"%s\", \"%s\");'%(dset,
                            dset.block_name, dtype, dset)
        temp = WriteString(temp)
        out = [WriteString("{"), halo_p, halo_m, sizes, base, value, temp, WriteString("}")]
        return out


    def add_block_name_to_kernel_sb(self, kernel):
        kernel.block_name = "block"
        return

    def generate_OPSC_dependants_sb(self, algorithm):
        def _generate(components):
            for component1 in components:
                if isinstance(component1, Loop):
                    return _generate(component1.components)
                elif isinstance(component1, Kernel):
                    self.Rational_constants = self.Rational_constants.union(component1.Rational_constants)
                    self.constants = self.constants.union(component1.constants).union(component1.IndexedConstants)

        code = algorithm.prg.opsc_code
        #print '\n'.join(code)
        _generate(algorithm.prg.components)
        return