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
import sympy.functions.special.tensor_functions as tf
import re
import inspect

import logging
LOG = logging.getLogger(__name__)


# Get Sympy Tensor Functions
SYMPY_FUNCTIONS = [str(m[0]) for m in inspect.getmembers(tf, inspect.isclass) if m[1].__module__ == tf.__name__]
LOCAL_FUNCTIONS = []


class Conservative(Derivative):
    LOCAL_FUNCTIONS.append('Conservative')


class Der(Derivative):
    LOCAL_FUNCTIONS.append('Der')
    def get_nested(self):
        """ Return any local classes nested within this derivative. """
        local_class, sympy_class = get_nested_classes(self)
        return local_class


class EinsteinTerm(Symbol):

    """ Represents any symbol in the equation as a SymPy Symbol object which in turn represents an Einstein term.
    This could be e.g. tau_i_j, but can also be e.g. u, rho.
    In other words, all symbols in the equation are Einstein terms, but they can have zero or more indices. """

    def __new__(self, symbol, **assumptions):
        self._sanitize(assumptions, self) # Remove any 'None's, etc.
        self.name = str(symbol)
        
        # Make this into a new SymPy Symbol object.
        self = Symbol.__xnew__(self, self.name, **assumptions)
        
        # Is this a term that is constant in space and time (i.e. doesn't have any indices).
        self.is_constant =  False

        # Extract the indices, which are always preceeded by an underscore.
        indices = self.name.split('_')[1:]
        self.indices = ['_' + x for x in indices]
        return self

    def get_indices(self):
        return self.indices

    def has_index(self, index):
        """ Check to see whether a given index exists in the Einstein variable. """
        found = False
        for i in self.indices:
            if i == index:
                found = True
                break
        return found

    def apply_index(self, index, value):
        """ Replace the index with a particular value (e.g. replace i with 0 in x_i), and return the updated EinsteinTerm object. """
        updated_symbol = str(self).replace(index, str(value))
        if self.is_constant:
            new = EinsteinTerm(updated_symbol)
            new.is_constant = True
            return new
        else:
            return EinsteinTerm(updated_symbol)
        return

def get_nested_classes(function):
    """ Match up the SymPy/local classes with the arguments of the function provided. """
    local_classes = []
    sympy_classes = []
    for arg in function.args:
        # Break down the arguments into individual Function-type objects.
        for atom in arg.atoms(Function):
            try:
                loc = LOCAL_FUNCTIONS.index('%s' % type(atom))
                local_classes.append(atom)
            except ValueError:
                # If the Function atom is not in our list of local/AutoFD function classes, then check to see whether it is a SymPy class instead.
                try:
                    loc = SYMPY_FUNCTIONS.index('%s' % type(atom))
                    sympy_classes.append(atom)
                except:
                    raise ValueError("The function provided (%s) is not a local nor a SymPy function type." % atom)

    return local_classes, sympy_classes


class EinsteinExpansion(object):

    """ Expand an Einstein variable with respect to its indices. """

    def __init__(self, expression, ndim):
        self.expression = expression
        self.ndim = ndim
        self.expanded = []
        
        self.is_equality = isinstance(expression, Equality)
        
        # If we have an equation, then expand the RHS only
        if self.is_equality:
            # Get LHS indices
            lhs_indices = self.get_einstein_indices(expression.lhs)
            # Get RHS indices
            rhs_indices = self.get_einstein_indices(expression.rhs)
            # Any non-common indices
            indices = rhs_indices.difference(lhs_indices)
            part_to_expand = expression.rhs
        else:
            indices = self.get_einstein_indices()
            part_to_expand = expression
        
        # Expand indices
        for index in indices:
            LOG.debug('Expanding with respect to index %s' % index)
            part_to_expand = self.apply(part_to_expand, index)
            part_to_expand = self.replace_sympy_functions(part_to_expand)

        # Apply the LHS if equality
        if self.is_equality:
            temp = Equality(expression.lhs, part_to_expand)
            if lhs_indices:
                vector_expansion = [temp for dim in range(self.ndim)]
                print(vector_expansion)
                for dim in range(self.ndim):
                    for index in lhs_indices:
                        for term in vector_expansion[dim].atoms(EinsteinTerm):
                            if term.has_index(index):
                                new = term.apply_index(index,dim)
                                vector_expansion[dim] = vector_expansion[dim].xreplace({term:new})
                    vector_expansion[dim] = self.replace_sympy_functions(vector_expansion[dim])
                self.expanded  = vector_expansion
                pprint(self.expanded)
            else:
                self.expanded = Equality(expression.lhs, part_to_expand)
        else:
            self.expanded = part_to_expand

    def replace_sympy_functions(self, part_to_apply):
        """ If the term(s) contain SymPy functions (e.g. Derivative), then replace them here. """
        local, sympy = get_nested_classes(part_to_apply)
        for function in sympy:
            index_exist = self.get_einstein_indices(function)
            if not index_exist:
                arg = [int(str(arg)) for arg in function.args]
                value = (type(function)(*arg))
                part_to_apply = part_to_apply.xreplace({function:value})
        return part_to_apply

    def has_nested_derivative(self, function):
        """ Return True if the function's arguments contains a Derivative, otherwise return False. """
        nested = False
        for arg in function.args:
            if arg.atoms(Derivative):
                nested = True
        return nested
    
    def check_index_nested(self, expression, index):
        """ This takes the expression as an input and returns the expression if all the terms in the expression
        has the Einstein index (index). Otherwise it returns a list of terms that has the Einstein Index. """

        has_index = []
        def _check_index(expression,index):
            arg_index = [False for arg in expression.args]
            for number,arg in enumerate(expression.args):
                if(any(term.has_index(index) for term in arg.atoms(EinsteinTerm))):
                    arg_index[number] = True
            nested = self.has_nested_derivative(expression)
            
            if all(arg_index):
                has_index.append(expression)
            elif not nested:
                if any(arg_index):
                    has_index.append(expression)
            elif any(arg_index[1:]):
                has_index.append(expression)
            elif nested:
                for number,arg in enumerate(expression.args):
                    if arg_index[number]:
                        if arg.is_Mul:
                            new_term =  self.check_index_MUL(arg,index)
                            
                            has_index.append(new_term)
                        elif arg.is_Add:
                            new_term = self.check_index_ADD(arg,index)
                            has_index.append(new_term)
                            
                        else:
                            local, sympy = get_nested_classes(expression)
                            for func in local:
                                _check_index(func,index)
                            for func in sympy:
                                _check_index(func,index)
            return
        _check_index(expression,index)
        if len(has_index) == 1 and has_index[0] == expression:
            return expression
        else:
            return flatten(has_index)
        return
    def check_index_ADD(self, add_term, index):
        terms = add_term.as_ordered_terms()
        is_indexed_add = [False for arg in terms]
        indexed_terms = []
        
        for number,term in enumerate(terms):
            if(any(term.has_index(index) for term in term.atoms(EinsteinTerm))):
                indexed_terms += [term]
        return indexed_terms

    def check_index_MUL(self,multiplication_term,index):
        terms = multiplication_term.as_ordered_factors()
        is_indexed_mul = [False for arg in terms]
        for number,term in enumerate(terms):
            if(any(term.has_index(index) for term in term.atoms(EinsteinTerm))):
                is_indexed_mul[number] = True
        out_term = multiplication_term
        for number,term in enumerate(terms):
            if isinstance (term,Symbol) and not term.has_index(index):
                out_term = out_term/term
            elif not is_indexed_mul[number]:
                out_term = out_term/term
        return out_term

    def get_einstein_indices(self, function=None):
        """ Return all the indices of a function that are Einstein indices. """
        
        # If no particular function is provided, then use the full expression provided when the expansion class was first set up.
        if function is None:
            expression = self.expression
        else:
            expression = function

        # Get all the atoms in the expression that are Einstein variables, and then return their indices.
        einstein_indices = []
        for atom in expression.atoms(EinsteinTerm):
            einstein_indices += atom.get_indices()

        return set(einstein_indices)

    def find_terms(self, terms, index):
        """ Return all terms with a given index. """

        final_terms = {}
        def _find_terms(term,index):
                if isinstance(term, Basic):
                    if term.is_Add:
                        args = self.check_index_ADD(term,index)
                        for arg in args:
                            _find_terms(arg,index)
                    elif term.is_Mul:
                        new_term =  self.check_index_MUL(term, index)
                        if term == new_term:
                            has_terms.append(term)
                        else:
                            _find_terms(new_term, index)
                    elif isinstance(term, Derivative):
                        new_function = self.check_index_nested(term,index)
                        
                        if term == new_function:
                            print("TERM EQUALS FUNCTION")
                            pprint(term); pprint(new_function)
                            has_terms.append(term)
                        else:
                            for f in new_function:
                                _find_terms(f, index)
                    elif isinstance(term, Function):
                        if any(term.has_index(index) for term in term.atoms(EinsteinTerm)):
                            has_terms.append(term)
                    elif isinstance(term, EinsteinTerm):
                        if term.has_index(index):
                            has_terms.append(term)
                    else:
                        # TODO raise Value Error
                        print('NOT ABLE TO Process', term, type(term.func))
                return
        for term, repr in terms:
            print("IN FIND TERMS", term)
            has_terms = []
            _find_terms(term,index)
            final_terms[term] = flatten(has_terms)
        return final_terms

    def expand_indices_terms(self, terms, index, dimensions):
        """ This expands the terms given in the list of terms

        :returns: a dictionary with each term as a key.
        :rtype: dict
        """

        expanded_terms = {}
        for term in terms:
            new_term = [term for dim in range(0,dimensions)]
            for at in term.atoms(EinsteinTerm):
                for dim in range(0,dimensions):
                    new = at.apply_index(index,dim)
                    new_term[dim] = new_term[dim].xreplace({at:new})
            expanded_terms[term] = new_term
        return expanded_terms

    def apply(self, part, index):
        ndim = self.ndim
        terms, gene = part.as_terms()
        final_terms = self.find_terms(terms, index)
        for term, repr in terms:
            if final_terms[term]:
                for new_term in final_terms[term]:
                    t = self.expand_indices_terms([new_term], index, self.ndim)
                    for key in t.keys():
                        val = sum(t[key])
                        part = part.xreplace({key:val})
                    print('Expanded term', term,t)
        return part

def apply_formulations(expression):
    '''
    Apply the formulations 
    # TODO Move applying Derivative for conservative or Der formulations to their respective
    classes and call the class from here
    
    '''
    temp_expression = expression
    der_dire = set()
    for atom in temp_expression.atoms(Derivative):
        for arg in atom.args[1:]:
            der_dire.add(arg)
    function_vars = set()
    for term in temp_expression.atoms(EinsteinTerm):
        function_vars.add(term)
    function_vars = function_vars.difference(der_dire)
    for term in function_vars:
        if not term.is_constant :
            temp_expression = temp_expression.replace(term, Function('%s'%term)(*der_dire))
    for atom in temp_expression.atoms(Der):
        new_atom = atom.doit()
        temp_expression = temp_expression.xreplace({atom:new_atom})
    temp_expression = temp_expression
    pprint(temp_expression)
    return temp_expression

class Equation(object):

    """ Describes an equation we want to solve. """

    def __init__(self, expression, problem):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :returns: None
        """
        local_dict = {'Der':Der, 'Symbol':EinsteinTerm, 'Conservative':Conservative} # TODO: automate from local classes

        self.original = expression

        # Perform substitutions, if any.
        if problem.substitutions:
            for sub in problem.substitutions:
                temp = parse_expr(sub, local_dict, evaluate = False)
                self.original = self.original.replace(str(temp.lhs), str(temp.rhs))

        # Parse the equation.
        self.parsed = parse_expr(self.original, local_dict, evaluate = False)
        self.expanded = []
        
        # Update the Einstein Variables in the expression that are constants to is_const = True
        for term in self.parsed.atoms(EinsteinTerm):
            if any(constant == str(term) for constant in problem.constants):
                term.is_constant = True

        # Expand Einstein terms/indices
        expansion = EinsteinExpansion(self.parsed, problem.ndim)
        if isinstance(expansion.expanded, list):
            for equation in expansion.expanded:
                self.expanded.append(apply_formulations(equation))
        else:
            self.expanded.append(apply_formulations(expansion.expanded))
        
        return

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
