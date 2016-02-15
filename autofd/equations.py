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
Sympy_tensor_function = [str(m[0]) for m in inspect.getmembers(tf, inspect.isclass) if m[1].__module__ == tf.__name__]
my_local_function = []


class EinstienVariable(Symbol):
    is_EinstienVariable = True
    is_constant = False
    def __new__(self, symbol,  **assumptions):
        self._sanitize(assumptions, self)
        name = str(symbol)
        self = Symbol.__xnew__(self, name,  **assumptions)
        indices = name.split('_')[1:]
        self.is_constant =  False
        self.is_conser = False
        self.indices = ['_' + x for x in indices]
        return self
    def get_indices(self):
        return self.indices
    def check_index(self,index):
        is_index = False
        for ind in self.indices:
            if ind == index:
                is_index = True
            else:
                pass
        return is_index
    def apply_index(self,index,value):
        new = str(self).replace(index,str(value))
        if self.is_constant:
            newvar = EinstienVariable(new)
            newvar.is_constant = True
            return newvar
        else:
            return EinstienVariable(new)
        return

def get_tensor_class(function):
    '''
    Get the classes of the function arguments
    '''
    local_class = []
    sympy_class = []
    for arg in function.args:
        for at in arg.atoms(Function):
            try:
                loc = my_local_function.index('%s'%type(at))
                local_class += [at]
                #print("LocalClass", str(type(at)))
            except:
                try:
                    loc = Sympy_tensor_function.index('%s'%type(at))
                    #print("Sympy", str(type(at)))
                    sympy_class += [at]
                except:
                    raise ValueError("Nothing")
    return local_class, sympy_class
def is_nested(function):
    nested = False
    for arg in function.args:
        if arg.atoms(Derivative):
            nested = True
    return nested

class Conservative(Derivative):
    my_local_function.append('Conservative')
    
class Der(Derivative):
    my_local_function.append('Der')
    def get_nested(self):
        '''
        Return if any local classes (Not required )
        '''
        local_class, sympy_class = get_tensor_class(self)
        
        return local_class

class EinsteinExpansion(object):
    '''
    This class apply the Einstein Expansion of an expression
    '''

    def __init__(self,expression, ndim):
        self.expression = expression
        self.ndim = ndim
        self.expanded = []
        # temporary equation
        temporary_expression = expression
        self.Equality = isinstance(temporary_expression, Equality)
        print(self.Equality)
        # get LHS indices
        if self.Equality:
            lhs_indices = self.get_Einstein_indices(temporary_expression.lhs)
            # Get RHS indices
            rhs_indices = self.get_Einstein_indices(temporary_expression.rhs)
            # Diff
            indices = rhs_indices.difference(lhs_indices)
        else:
            indices = self.get_Einstein_indices()
        if self.Equality:
            part_to_expand = temporary_expression.rhs
        else:
            part_to_expand = temporary_expression
        # expand indices
        for index in indices:
            print('EXPNADING IN index', index)
            part_to_expand = self.apply_Einstien_expansion(part_to_expand,index)
            part_to_expand = self.apply_sympy_function(part_to_expand)
        # Apply the LHS if equality
        if self.Equality:
            temp = Equality(temporary_expression.lhs,part_to_expand)
            if lhs_indices:
                vector_expansion = [temp for dim in range(self.ndim)]
                print(vector_expansion)
                for dim in range(self.ndim):
                    for index in lhs_indices:
                        for Ev in vector_expansion[dim].atoms(EinstienVariable):
                            if Ev.check_index(index):
                                new = Ev.apply_index(index,dim)
                                vector_expansion[dim] = vector_expansion[dim].xreplace({Ev:new})
                    vector_expansion[dim] = self.apply_sympy_function(vector_expansion[dim])
                self.expanded  = vector_expansion
                pprint(self.expanded)
            else:
                self.expanded = Equality(temporary_expression.lhs,part_to_expand)
        else:
            self.expanded = part_to_expand

    def apply_sympy_function(self, part_to_apply):
        local, sympy = get_tensor_class(part_to_apply)
        for fn in sympy:
            index_exist = self.get_Einstein_indices(fn)
            if not index_exist:
                arg = [int(str(arg)) for arg in fn.args]
                value = (type(fn)(*arg))
                part_to_apply = part_to_apply.xreplace({fn:value})
        return part_to_apply

    def check_index_nested(self, expression, index):
        '''
        This takes the Function as an input and returns the expression if all the terms in the expression
        has the Einstein index (index).

        Else it returns a list of terms that has the Einstein Index
        '''
        has_index = []
        def _check_index(expression,index):
            arg_index = [False for arg in expression.args]
            for number,arg in enumerate(expression.args):
                if(any(Ev.check_index(index) for Ev in arg.atoms(EinstienVariable))):
                    arg_index[number] = True
            nested = is_nested(expression)
            
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
                            newt =  self.check_index_MUL(arg,index)
                            
                            has_index.append(newt)
                        elif arg.is_Add:
                            newt = self.check_index_ADD(arg,index)
                            has_index.append(newt)
                            
                        else:
                            local, sympy = get_tensor_class(expression)
                            for func in local:
                                _check_index(func,index)
                            for func in sympy:
                                _check_index(func,index)
            return
        _check_index(expression,index)
        if len(has_index) == 1 and has_index[0] == expression:
            #print()
            return expression
        else:
            return flatten(has_index)
        return
    def check_index_ADD(self, add_term, index):
        terms = add_term.as_ordered_terms()
        is_indexed_add = [False for arg in terms]
        indexed_terms = []
        
        for number,term in enumerate(terms):
            if(any(Ev.check_index(index) for Ev in term.atoms(EinstienVariable))):
                indexed_terms += [term]
        #print("IN CHECK Add", add_term, terms)
        #print("AND THE TERMS ARE", indexed_terms)
        return indexed_terms

    def check_index_MUL(self,multiplication_term,index):
        terms = multiplication_term.as_ordered_factors()
        is_indexed_mul = [False for arg in terms]
        for number,term in enumerate(terms):
            if(any(Ev.check_index(index) for Ev in term.atoms(EinstienVariable))):
                is_indexed_mul[number] = True
        out_term = multiplication_term
        for number,term in enumerate(terms):
            if isinstance (term,Symbol) and not term.check_index(index):
                #print("Not", Symbol)
                out_term = out_term/term
            elif not is_indexed_mul[number]:
                out_term = out_term/term
        #print("MULS INDEXED", multiplication_term, out_term)
        return out_term

    def get_Einstein_indices(self, function=None):
        if function == None:
            equation = self.expression
        else:
            equation = function
        Einstein_indices = []
        for atom in equation.atoms(EinstienVariable):
            Einstein_indices += atom.get_indices()
        return set(Einstein_indices)

    def find_terms(self,terms,index):
        final_terms = {}
        def _find_terms(term,index):
                if isinstance(term, Basic):
                    if term.is_Add:
                        args = self.check_index_ADD(term,index)
                        for arg in args:
                            _find_terms(arg,index)
                    elif term.is_Mul:
                        newt =  self.check_index_MUL(term,index)
                        if term == newt:
                            has_terms.append(term)
                        else:
                            _find_terms(newt,index)
                    elif isinstance(term, Derivative):
                        #print(term)
                        newF = self.check_index_nested(term,index)
                        
                        if term == newF:
                            print("TERM EQUALS FUNCTION")
                            pprint(term); pprint(newF)
                            has_terms.append(term)
                        else:
                            for f in newF:
                                _find_terms(f,index)
                    elif isinstance(term, Function):
                        if any(Ev.check_index(index) for Ev in term.atoms(EinstienVariable)):
                            has_terms.append(term)
                    elif isinstance(term, EinstienVariable):
                        if term.check_index(index):
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
        return  final_terms

    def expand_indices_terms(self,terms,index, dimensions):
        '''
        This expands the terms given in the list of terms

        returns a dictionary with each term as a key

        '''
        expanded_terms ={term:None for term in terms}
        for term in terms:
            new_term = [term for dim in range(0,dimensions)]
            for at in term.atoms(EinstienVariable):
                for dim in range(0,dimensions):
                    new = at.apply_index(index,dim)
                    new_term[dim] = new_term[dim].xreplace({at:new})
            expanded_terms[term] = new_term
        return expanded_terms

    def apply_Einstien_expansion(self,part,index):
        ndim = self.ndim
        terms, gene = part.as_terms()
        final_terms = self.find_terms(terms,index)
        #print("Input part", part)
        #print("FT",final_terms)
        for term,repr in terms:
            if final_terms[term]:
                for new_term in final_terms[term]:
                    t = self.expand_indices_terms([new_term],index, ndim)
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
    for ev in temp_expression.atoms(EinstienVariable):
        function_vars.add(ev)
    function_vars = function_vars.difference(der_dire)
    for ev in function_vars:
        if not ev.is_constant :
            temp_expression = temp_expression.replace(ev, Function('%s'%ev)(*der_dire))
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
        local_dict = {'Der':Der,'Symbol':EinstienVariable,'Conservative':Conservative} # automate from local classes

        self.original = expression
        #self.ndim = problem.ndim

        # Perform substitutions, if any.
        if problem.substitutions:
            for sub in problem.substitutions:
                temp = parse_expr(sub,local_dict,evaluate = False)
                self.original = self.original.replace(str(temp.lhs), str(temp.rhs))

        # Parse the equation.
        self.parsed = parse_expr(self.original,local_dict,evaluate = False)
        self.expanded = []
        
        # Update the Einstein Variables in the expression that are constants to is_const = True
        for Ev in self.parsed.atoms(EinstienVariable):
            if any(const == str(Ev) for const in problem.constants):
                Ev.is_constant = True
                print(Ev)
        # Expand Einstein Indices
        EE_expression = EinsteinExpansion(self.parsed,problem.ndim)
        if isinstance(EE_expression.expanded, list):
            for eq in EE_expression.expanded:
                self.expanded.append(apply_formulations(eq))
        else:
            self.expanded.append(apply_formulations(EE_expression.expanded))
        
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
