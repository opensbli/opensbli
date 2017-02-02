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
class OPSC(object):

    def __init__(self, algorithm):
        """ Generating an OPSC code from the algorithm"""
        if not algorithm.MultiBlock:
            self.datasets_to_declare = []
            self.Rational_constants = set()
            self.constants = set()
            self.stencils_to_declare = set()
            self.generate_OPSC_dependants_sb(algorithm)
        return

    def add_block_name_to_kernel_sb(self, kernel):
        kernel.block_name = "block"
        return
    def kernel_datasets(self, kernel):
        lhs = list(kernel.lhs_datasets) +list(kernel.rhs_datasets)
        #rhs = kernel.rhs_datasets
        dsets = []
        for d in lhs:
            d1 = DeclareDataset(d, kernel.block_number)
            print d, d1
        return dsets


    def generate_OPSC_dependants_sb(self, algorithm):
        def _generate(components):
            for component1 in components:
                if isinstance(component1, Loop):
                    return _generate(component1.components)
                elif isinstance(component1, Kernel):
                    self.add_block_name_to_kernel_sb(component1)
                    #self.kernel_datasets(component1)
                    self.datasets_to_declare += list(component1.lhs_datasets) + list(component1.rhs_datasets)
                    self.Rational_constants = self.Rational_constants.union(component1.Rational_constants)
                    self.constants = self.constants.union(component1.constants).union(component1.IndexedConstants)
                    # WARNING need to implement this
                    stens = component1.get_stencils
                    for key, value in stens.iteritems():
                        self.stencils_to_declare.add(tuple(value))
        code = algorithm.prg.opsc_code
        #print '\n'.join(code)
        _generate(algorithm.prg.components)
        return