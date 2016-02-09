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
from sympy.parsing.sympy_parser import parse_expr
import re

import logging
LOG = logging.getLogger(__name__)


class Equation(object):

    """ Describes an equation we want to solve. """

    def __init__(self, equation, system):
        self.original = equation

        # Perform substitutions, if any
        if system.substitutions:
            for sub in system.substitutions:
                temp = parse_expr(sub)
                self.original = self.original.replace(str(temp.lhs), str(temp.rhs))

        # Parse the equation
        self.parsed = parse_expr(self.original)
        self.expanded = []

        index = find_indices(self.parsed.lhs.atoms(Symbol))
        ndim = system.ndim
        if index:
            for dim in range(0, system.ndim):
                for ind in index:
                    lhs = parse_expr(str(self.parsed.lhs).replace(str(ind), str(dim)))
                    rhs = parse_expr(str(self.parsed.rhs).replace(str(ind), str(dim)))
                    e = self.parsed.replace(self.parsed.lhs, lhs).replace(self.parsed.rhs, rhs)
                    self.expanded = self.expanded + [e]
        else:
            self.expanded = self.expanded + [self.parsed]

        self.constants = list(Symbol(c) for c in system.constants)

        # Treat special operators such as conser and skew.
        for equation_number in range(len(self.expanded)):
            e = self.expanded[equation_number]
            indices = find_indices(e.rhs.atoms(Symbol))

            all_operators = list(e.atoms(Function('conser'))) + list(e.atoms(Function('Der'))) + list(e.atoms(Function('Skew')))

            if all_operators:
                operator_args = flatten(list(set(o.args[1:] for o in all_operators)))
            else:
                operator_args = []
            if indices:
                syms = list(flatten(list(e.atoms(Symbol).difference(set(self.constants + operator_args + list(Symbol(ind) for ind in indices))))))
                if operator_args:
                    temp = list(str(o) for o in operator_args)
                    fndef = ','.join(temp)
                    fns = []
                    for atom in self.expanded[equation_number].atoms(Function('conser')):
                        out = conser(atom)
                        self.expanded[equation_number] = self.expanded[equation_number].replace(atom, out)
                    for atom in self.expanded[equation_number].atoms(Function('Skew')):
                        out = skew(atom)
                        self.expanded[equation_number] = self.expanded[equation_number].replace(atom, out)
                    for sym in syms:
                        fn = parse_expr('%s(%s)' % (sym, fndef))
                        fns = fns + [fn]
                        self.expanded[equation_number] = self.expanded[equation_number].replace(sym, fn)

                    self.expanded[equation_number] = der(self.expanded[equation_number])
                    for no in range(len(syms)):
                        self.expanded[equation_number] = self.expanded[equation_number].subs(fns[no], syms[no])
                
                # Find the terms containing index i.
                for i in indices:
                    terms = find_terms(self.expanded[equation_number].rhs, [i])
                    for term in terms:
                        self.expanded[equation_number] = expand_ind(term, system.ndim, [i], self.expanded[equation_number])

        self.variables = []
        self.conser = []
        for equation_number, e in enumerate(self.expanded):
            e = self.expanded[equation_number]
            derivatives = list(e.lhs.atoms(Derivative))
            if len(derivatives) > 1:
                raise ValueError('More than one derivative in LHS')
            elif len(derivatives) == 0:
                pass
            else:
                self.conser.append(derivatives[0].args[0])
            all_derivatives = list(e.atoms(Derivative))
            if all_derivatives:
                derivative_args = flatten(list(set(d.args[1:] for d in all_derivatives)))
            else:
                derivative_args = []
            symvar = list(flatten(list(e.atoms(Symbol).difference(set(self.constants + derivative_args + self.conser)))))
            self.variables.append(symvar)

            self.const = list(flatten(list(e.atoms(Symbol).difference(set(self.variables[equation_number] + derivative_args + self.conser)))))
            self.ndim = system.ndim

        return

def find_terms(expr, indices):
    terms = []
    if any(str(expr).count(str(ind)) for ind in indices):
        if expr.is_Mul:
            facs = expr.as_two_terms()
            # for fac in facs:
            if all(any(str(fac).count(str(ind)) for ind in indices) for fac in facs):
                terms = expr
            else:
                if any(fac.is_Mul or fac.is_Add for fac in facs):
                    for fac in facs:
                        if any(str(fac).count(str(ind)) for ind in indices):
                            find_terms(fac, indices)
                else:
                    for fac in facs:
                        if any(str(fac).count(str(ind)) for ind in indices):
                            find_terms(fac, indices)
        elif expr.is_Add:
            facs = expr.as_two_terms()
            if all(any(str(fac).count(str(ind)) for ind in indices) for fac in facs):
                if any(fac.is_Mul or fac.is_Add for fac in facs):
                    for fac in facs:
                        find_terms(fac, indices)
                else:
                    terms = expr
            else:
                for fac in facs:
                    find_terms(fac, indices)

        else:
            terms = expr
    else:
        terms = []

    return terms

def conser(inp):
    out = parse_expr('Derivative(%s,%s)' % (inp.args[0], inp.args[1]))
    return out

def skew(inp):
    global indterm
    # This is to be done
    mult = inp.args[0].as_two_terms()
    pprint(mult)
    te1 = 'Derivative(%s,%s)' % (inp.args[0], inp.args[1])
    te2 = '(%s) * Derivative(%s,%s)' % (mult[0], mult[1], inp.args[1])
    te3 = '(%s) * Derivative(%s,%s)' % (mult[1], mult[0], inp.args[1])
    out = parse_expr('0.5*(%s + %s + %s)' % (te1, te2, te3))
    pprint(out)
    indice = find_indices(out.atoms(Symbol))
    for ind in indice:
        indterm = []
        index = [ind]
        find_terms(out, index)
        for term in indterm:
            out = expand_ind(term, ndim, index, out)
    return out

def der(inp):
    lhs = parse_expr(str(inp.lhs).replace('Der(', 'diff('))
    rhs = parse_expr(str(inp.rhs).replace('Der(', 'diff('))
    out = Eq(lhs, rhs)
    return out
            
class ExplicitFilter(object):

    """ Explicit filter implemented from the NASA PAPER. """

    def __init__(self):
        self.coeffs = {}
        self.coeffs[2] = [-1, 2, -1]
        self.coeffs[4] = [-1, 4, -6, 4, -1]

    def filter_equations(self, alg):
        """ Returns the filter equations. """
        filtereqs = []
        for dire, fil in enumerate(alg.expfilter):
            if fil:
                ooa = alg.spatial_order
                alphaD = (-1)**(int(ooa/2) + 1) * (2) ** (-ooa)
                if self.coeffs.get(ooa):
                    muls = self.coeffs.get(ooa)
                    temp = list(con.base for con in inp.conser)
                    out = [con for con in inp.conser]
                    # out = 0
                    direction = Symbol('x%d' % dire)
                    for num, ind in enumerate(inp.fd_points):
                        repl = '%s' % str(direction + ind)
                        index = str(inp.varform).replace(str(direction), repl)
                        index = Idx(index, Symbol('nblock', integer=True))
                        for nv in range(len(temp)):
                            varib = temp[nv]
                            out[nv] += alphaD * muls[num] * varib[index]
                    for nv in range(len(temp)):
                        out[nv] = Eq(savedic.get(inp.conser[nv]), out[nv])
                    filtereqs.append(out)
            else:
                raise NotImplementedError('Implement spatial filtering coefficients for order of accuracy %d' % ooa)
        return filtereqs


def find_indices(equations):
    """ Finds the variable indices in all the input equations and returns a list of the set of indices.

    :arg equations: a string, symbolic equation, a list of strings, or a list of symbolic equations.
    :returns: a list of the set of indices.
    :rtype: list
    """

    # All indices after after an underscore, so match this using a regular expression.
    index = []
    if isinstance(equations, list) or isinstance(equations, set):
        for e in equations:
            e = str(e)
            index = index + list(set(re.compile(r'\_\S', re.I).findall(e)))
    else:
        e = str(equations)
        index = index + list(set(re.compile(r'\_\S', re.I).findall(e)))
    index = list(set(index))
    return index


def expand_ind(term, ndim, index, equation):
    """ Utility function to expand a given equations """

    fin = ''
    if term.is_Mul:
        out = str(term)
        # index = tensor_indices(self)
        for ind in index:
            for dim in range(0, int(ndim)):
                fin = fin + '+' + str(out).replace(str(ind), str(dim))
            out = fin
            fin = ''
    elif term.is_Add:
        fin1 = ''
        for te in term.as_ordered_terms():
            fin = ''
            # pprint(te)
            out = str(te)
            # index = tensor_indices(te)
            for ind in index:
                for dim in range(0, int(ndim)):
                    fin = fin + '+' + str(out).replace(str(ind), str(dim))
                out = fin
                fin = ''
            fin1 = fin1 + '+' + out
        out = fin1
    else:
        out = str(term)
        # index = tensor_indices(self)
        for ind in index:
            for dim in range(0, int(ndim)):
                fin = fin + '+' + str(out).replace(str(ind), str(dim))
            out = fin
            fin = ''
    # pprint(out)
    out = parse_expr(out)
    # pprint(out)
    equation = equation.subs(term, out)
    return equation


def variable_count(variable, equations):
    """ Return the number of input equations containing a particular variable.
    
    :arg variable: the variable under consideration.
    :arg equations: the equations to search.
    :returns: the number of equations containing the variable.
    :rtype: int
    """

    count = 0
    for e in equations:
        count = count + e.count(variable)
    return count


def equations_to_dict(equations):
    """ Get the LHS and RHS of each equation, and return them in dictionary form.
    
    :arg equations: the equations to consider.
    :returns: a dictionary of (LHS, RHS) pairs.
    :rtype: dict
    """

    lhs = list(e.lhs for e in equations)
    rhs = list(e.rhs for e in equations)
    d = dict(zip(lhs, rhs))
    return d
