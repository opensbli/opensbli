#!/usr/bin/env python

from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import re

import logging
LOG = logging.getLogger(__name__)


class Equation(object):

    """ Describes an equation we want to solve. """

    global indterm

    def __init__(self, eq, system):
        self.inpeq = eq

        # Perform substitutions if any
        if system.substitutions:
            for form in system.substitutions:
                temp = parse_expr(form)
                self.inpeq = self.inpeq.replace(str(temp.lhs), str(temp.rhs))

        # Parse the equation
        self.parsed = parse_expr(self.inpeq)
        self.expandedeq = []

        index = tensor_indices(self.parsed.lhs.atoms(Symbol))
        ndim = system.ndim
        if index:
            for dim in range(0, system.ndim):
                for ind in index:
                    lhs = parse_expr(str(self.parsed.lhs).replace(str(ind), str(dim)))
                    rhs = parse_expr(str(self.parsed.rhs).replace(str(ind), str(dim)))
                    eqn = self.parsed.replace(self.parsed.lhs, lhs).replace(self.parsed.rhs, rhs)
                    self.expandedeq = self.expandedeq + [eqn]
        else:
            self.expandedeq = self.expandedeq + [self.parsed]
        self.constants = list(Symbol(con) for con in system.constants)

        def fin_terms(expr, indices):
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
                                    fin_terms(fac, indices)
                        else:
                            for fac in facs:
                                if any(str(fac).count(str(ind)) for ind in indices):
                                    fin_terms(fac, indices)
                elif expr.is_Add:
                    facs = expr.as_two_terms()
                    if all(any(str(fac).count(str(ind)) for ind in indices) for fac in facs):
                        if any(fac.is_Mul or fac.is_Add for fac in facs):
                            for fac in facs:
                                fin_terms(fac, indices)
                        else:
                            terms = expr
                    else:
                        for fac in facs:
                            fin_terms(fac, indices)

                else:
                    terms = expr
            else:
                terms = []
            if terms:
                indterm.append(terms)
            # pprint(terms)
            return

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
            indice = tensor_indices(out.atoms(Symbol))
            for ind in indice:
                indterm = []
                index = [ind]
                fin_terms(out, index)
                for term in indterm:
                    out = expand_ind(term, ndim, index, out)
            return out

        def der(inp):
            lhs = parse_expr(str(inp.lhs).replace('Der(', 'diff('))
            rhs = parse_expr(str(inp.rhs).replace('Der(', 'diff('))
            out = Eq(lhs, rhs)
            return out

        global indterm
        indterm = []

        for eqno in range(len(self.expandedeq)):
            eq = self.expandedeq[eqno]
            indices = tensor_indices(eq.rhs.atoms(Symbol))
            allfn = list(eq.atoms(Function('conser'))) + list(eq.atoms(Function('Der'))) + list(eq.atoms(Function('Skew')))
            if allfn:
                fn_args = flatten(list(set(der.args[1:] for der in allfn)))
            else:
                fn_args = []
            if indices:
                syms = list(flatten(list(eq.atoms(Symbol).difference(set(self.constants + fn_args + list(Symbol(ind) for ind in indices))))))
                if fn_args:
                    temp = list(str(fn) for fn in fn_args)
                    fndef = ','.join(temp)
                    fns = []
                    for atom in self.expandedeq[eqno].atoms(Function('conser')):
                        out = conser(atom)
                        self.expandedeq[eqno] = self.expandedeq[eqno].replace(atom, out)
                    for atom in self.expandedeq[eqno].atoms(Function('Skew')):
                        out = skew(atom)
                        self.expandedeq[eqno] = self.expandedeq[eqno].replace(atom, out)
                    for sym in syms:
                        fn = parse_expr('%s(%s)' % (sym, fndef))
                        fns = fns + [fn]
                        self.expandedeq[eqno] = self.expandedeq[eqno].replace(sym, fn)

                    self.expandedeq[eqno] = der(self.expandedeq[eqno])
                    for no in range(len(syms)):
                        self.expandedeq[eqno] = self.expandedeq[eqno].subs(fns[no], syms[no])
                for ind in indices:
                    indterm = []
                    index = [ind]
                    fin_terms(self.expandedeq[eqno].rhs, index)
                    # pprint(ind)
                    # pprint(indterm)
                    for term in indterm:
                        self.expandedeq[eqno] = expand_ind(term, system.ndim, index, self.expandedeq[eqno])

        self.variables = []
        self.conser = []
        for eqno, eq in enumerate(self.expandedeq):
            eq = self.expandedeq[eqno]
            ders = list(eq.lhs.atoms(Derivative))
            if len(ders) > 1:
                raise ValueError('More than one derivative in LHS')
            elif len(ders) == 0:
                pass
            else:
                self.conser.append(ders[0].args[0])
            allfn = list(eq.atoms(Derivative))
            if allfn:
                fn_args = flatten(list(set(der.args[1:] for der in allfn)))
            else:
                fn_args = []
            symvar = list(flatten(list(eq.atoms(Symbol).difference(set(self.constants + fn_args + self.conser)))))
            self.variables.append(symvar)

            self.const = list(flatten(list(eq.atoms(Symbol).difference(set(self.variables[eqno] + fn_args + self.conser)))))
            # pprint(eq)
            self.ndim = system.ndim

        return

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


def tensor_indices(equations):
    """ Finds the tensor indices in all the input equations and returns a list of the set of indices.

    :arg equations: a string, symbolic equations, a list of strings, or a list of symbolic equations.
    :returns: a list of the set of indices.
    :rtype: list
    """

    index = []
    if isinstance(equations, list):
        for eq in equations:
            eq = str(eq)
            index = index + list(set(re.compile(r'\_\S', re.I).findall(eq)))
    else:
        eq = str(equations)
        index = index + list(set(re.compile(r'\_\S', re.I).findall(eq)))
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
    LOG.debug(d)
    return d
