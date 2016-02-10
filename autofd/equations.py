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
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :returns: None
        """

        self.original = equation
        self.ndim = system.ndim

        # Perform substitutions, if any.
        if system.substitutions:
            for sub in system.substitutions:
                temp = parse_expr(sub)
                self.original = self.original.replace(str(temp.lhs), str(temp.rhs))

        # Parse the equation.
        self.parsed = parse_expr(self.original)
        self.expanded = []

        indices = find_indices(self.parsed.lhs.atoms(Symbol))
        if indices:
            for dim in range(0, system.ndim):
                for i in indices:
                    lhs = parse_expr(str(self.parsed.lhs).replace(str(i), str(dim)))
                    rhs = parse_expr(str(self.parsed.rhs).replace(str(i), str(dim)))
                    e = self.parsed.replace(self.parsed.lhs, lhs).replace(self.parsed.rhs, rhs)
                    self.expanded = self.expanded + [e]
        else:
            self.expanded = self.expanded + [self.parsed]

        # Get all constants and convert them to Symbols.
        self.constants = list(Symbol(c) for c in system.constants)

        # Treat special operators such as conser and skew.
        for equation_number in range(len(self.expanded)):
            e = self.expanded[equation_number]
            
            # Get a list of all function calls in the equation under consideration.
            all_operators = list(e.atoms(Function('conser'))) + list(e.atoms(Function('Der'))) + list(e.atoms(Function('Skew')))
            if all_operators:
                # Get the list of arguments for all the function calls.
                operator_args = flatten(list(set(o.args[1:] for o in all_operators)))
            else:
                operator_args = []

            # Find the index (e.g. i and j) of each symbol (e.g. x and y).
            indices = find_indices(e.rhs.atoms(Symbol))
            if indices:
                symbols = list(flatten(list(e.atoms(Symbol).difference(set(self.constants + operator_args + list(Symbol(index) for index in indices))))))
                if operator_args:
                    temp = list(str(o) for o in operator_args)
                    args = ','.join(temp)
                    functions = []
                    # conser
                    for atom in self.expanded[equation_number].atoms(Function('conser')):
                        out = conser(atom)
                        self.expanded[equation_number] = self.expanded[equation_number].replace(atom, out)
                    # Skew
                    for atom in self.expanded[equation_number].atoms(Function('Skew')):
                        out = skew(atom)
                        self.expanded[equation_number] = self.expanded[equation_number].replace(atom, out)
                    # All other function calls
                    for s in symbols:
                        f = parse_expr('%s(%s)' % (s, args))
                        functions += [f]
                        self.expanded[equation_number] = self.expanded[equation_number].replace(s, f)

                    self.expanded[equation_number] = der(self.expanded[equation_number])
                    for i in range(len(symbols)):
                        self.expanded[equation_number] = self.expanded[equation_number].subs(functions[i], symbols[i])
                
                # Find the terms containing index i.
                for i in indices:
                    terms = find_terms(self.expanded[equation_number].rhs, [i])
                    for term in terms:
                        self.expanded[equation_number] = expand_indices(term, system.ndim, [i], self.expanded[equation_number])

        # Expand derivatives
        self.variables = []
        self.conser = []
        for equation_number, e in enumerate(self.expanded):
            e = self.expanded[equation_number]
            # Derivatives on the LHS
            derivatives = list(e.lhs.atoms(Derivative))
            if len(derivatives) > 1:
                raise ValueError('More than one derivative in LHS')
            elif len(derivatives) == 0:
                pass
            else:
                self.conser.append(derivatives[0].args[0])
            
            # All derivatives
            all_derivatives = list(e.atoms(Derivative))
            if all_derivatives:
                derivative_args = flatten(list(set(d.args[1:] for d in all_derivatives)))
            else:
                derivative_args = []
            symvar = list(flatten(list(e.atoms(Symbol).difference(set(self.constants + derivative_args + self.conser)))))
            self.variables.append(symvar)

        return

def find_terms(expression, indices):
    """ Find all terms in a given expression with given indices.
    
    :arg str expression: the string expression containing the equation to consider
    :arg list indices: the list of indices (represented as strings) to look out for
    :returns: a list of terms with particular indices
    :rtype: list
    """

    terms = []
    if any(str(expression).count(str(index)) for index in indices):
        # Multiplication
        if expression.is_Mul:
            parts = expression.as_two_terms()
            if all(any(str(part).count(str(index)) for index in indices) for part in parts):
                terms = [expression]
            else:
                if any(part.is_Mul or part.is_Add for part in parts):
                    for part in parts:
                        if any(str(part).count(str(index)) for index in indices):
                            found = find_terms(part, indices)
                            if found:
                                terms += found
                else:
                    for part in parts:
                        if any(str(part).count(str(index)) for index in indices):
                            found = find_terms(part, indices)
                            if found:
                                terms += found
        # Addition
        elif expression.is_Add:
            parts = expression.as_two_terms()
            if all(any(str(part).count(str(index)) for index in indices) for part in parts):
                if any(part.is_Mul or part.is_Add for part in parts):
                    for part in parts:
                        found = find_terms(part, indices)
                        if found:
                            terms += found
                else:
                    terms = [expression]
            else:
                for part in parts:
                    found = find_terms(part, indices)
                    if found:
                        terms += found

        else:
            terms = [expression]

    return terms

def conser(inp):
    out = parse_expr('Derivative(%s,%s)' % (inp.args[0], inp.args[1]))
    return out

def skew(inp):
    # FIXME: This is to be done
    mult = inp.args[0].as_two_terms()
    pprint(mult)
    te1 = 'Derivative(%s,%s)' % (inp.args[0], inp.args[1])
    te2 = '(%s) * Derivative(%s,%s)' % (mult[0], mult[1], inp.args[1])
    te3 = '(%s) * Derivative(%s,%s)' % (mult[1], mult[0], inp.args[1])
    out = parse_expr('0.5*(%s + %s + %s)' % (te1, te2, te3))
    pprint(out)
    indices = find_indices(out.atoms(Symbol))
    for i in indices:
        found = find_terms(out, [i])
        for term in found:
            out = expand_indices(term, ndim, [i], out)
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


def expand_indices(term, ndim, indices, equation):
    """ Expand a given term with respect to its indices (e.g. expand u_i to u0, u1, u2 in 3D).
    Once finished, replace the expanded terms in the equation provided. """

    temp = ''
    if term.is_Mul:
        expanded = str(term)
        for i in indices:
            for dim in range(0, int(ndim)):
                temp = temp + '+' + str(expanded).replace(str(i), str(dim))
            expanded = temp
            temp = ''

    elif term.is_Add:
        temp1 = ''
        for t in term.as_ordered_terms():
            temp2 = ''
            expanded = str(t)
            for i in indices:
                for dim in range(0, int(ndim)):
                    temp2 = temp2 + '+' + str(expanded).replace(str(i), str(dim))
                expanded = temp2
                temp2 = ''
            temp1 = temp1 + '+' + expanded
        expanded = temp1

    else:
        expanded = str(term)
        for i in indices:
            for dim in range(0, int(ndim)):
                temp = temp + '+' + str(expanded).replace(str(i), str(dim))
            expanded = temp
            temp = ''

    expanded = parse_expr(expanded)
    equation = equation.subs(term, expanded)

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
