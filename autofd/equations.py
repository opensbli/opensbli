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
from sympy.core.assumptions import ManagedProperties
from sympy.parsing.sympy_parser import parse_expr
import sympy.functions.special.tensor_functions as tf
import re
import sys
from .array import MutableDenseNDimArray,  derive_by_array
import inspect
import sympy.core as core

import logging
LOG = logging.getLogger(__name__)

# Get Sympy Tensor Functions
SYMPY_FUNCTIONS = [str(m[0]) for m in inspect.getmembers(tf, inspect.isclass) if m[1].__module__ == tf.__name__]
LOCAL_FUNCTIONS = []


class Conservative(Function):
    LOCAL_FUNCTIONS.append('Conservative')
    is_commutative = False

class KD(Function):
    #LOCAL_FUNCTIONS.append('Der')
    is_commutative = False

class Der(Function):
    LOCAL_FUNCTIONS.append('Der')
    is_commutative = False
    
    def doit(self):
        return diff(*args)


class EinsteinTerm(Symbol):

    """ Represents any symbol in the equation as a SymPy Symbol object which in turn represents an Einstein term.
    This could be e.g. tau_i_j, but can also be e.g. u, rho.
    In other words, all symbols in the equation are Einstein terms, but they can have zero or more indices. """
    is_commutative = False
    def __new__(self, symbol, **assumptions):
        self._sanitize(assumptions, self) # Remove any 'None's, etc.
        self.name = str(symbol)

        # Make this into a new SymPy Symbol object.
        self = Symbol.__xnew__(self, self.name, **assumptions)

        # Is this a term that is constant in space and time (i.e. doesn't have any indices).
        self.is_constant = False

        # Extract the indices, which are always preceeded by an underscore.
        indices = self.name.split('_')[1:]
        self.indices = [Symbol(x, integer=True) for x in indices]
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

    def expand_index(self, index, value):
        """ Replace the index with a particular value (e.g. replace i with 0 in x_i), and return the updated EinsteinTerm object. """
        updated_symbol = str(self).replace(index, str(value))
        if self.is_constant:
            new = EinsteinTerm(updated_symbol)
            new.is_constant = True
            return new
        else:
            return EinsteinTerm(updated_symbol)



class EinsteinExpansion(object):

    """ Expand an Einstein variable with respect to its indices. """

    def __init__(self, expression, ndim):
        self.expression = expression
        self.expression_indexed = expression
        self.ndim = ndim
        self.expanded = []

        self.is_equality = isinstance(expression, Equality)
        fns = self.get_functions(expression)
        self.dictionary = {}
        self.indexed_object_no = 0
        self.indexedObject_name = 'Arr'
        for ind,fn in enumerate(fns):
            dictionary = {}
            print('\n');
            if len(fn.atoms(Function)) == 1:
                temp,dictionary = self.funtion_to_indexed(fn,dictionary)
                #print(temp)
                self.expression_indexed = self.expression_indexed.xreplace({fn:temp})
                print "SELF DICTIONARY BEFOR UPDATE", self.dictionary
                self.update_dictionary(dictionary)
                print '\n'
                print "Returned dictionary", dictionary
                print "SELF DICTIONARY", self.dictionary
            else:
                temp,dictionary = self.nested_function_to_indexed(fn,dictionary)
                self.expression_indexed = self.expression_indexed.xreplace({fn:temp})
                print "SELF DICTIONARY BEFOR UPDATE", self.dictionary
                self.update_dictionary(dictionary)
                print '\n'
                print "Returned dictionary", dictionary
                print "SELF DICTIONARY", self.dictionary
                #pprint(fn)
        print('\n')
        print "The input equation is :: ", self.expression
        print "The updated equation with indexed objects is :: ", self.expression_indexed

        print "The LHS Indices of the updated equation are :: ", get_indices(self.expression_indexed.lhs)

        print "The RHS Indices of the updated equation are :: ", get_indices(self.expression_indexed.rhs)

        print "Are LHS indices and rhs indices are same :: ", (get_indices(self.expression_indexed.lhs) == get_indices(self.expression_indexed.rhs))

        print('\n')
        print "The dictionary of evaluations for the functions is "
        for key,value in self.dictionary.iteritems():
            print key, ':', value
        print "Contraction structue is ", get_contraction_structure(self.expression_indexed.rhs)
        pprint(get_contraction_structure(self.expression_indexed.rhs))
        
        # now processing each arrays in the dictionary
        arrays = {}
        x = [symbols('x%d'%dim) for dim in range(ndim)]
        u = [Function('u%d'%dim)(*x) for dim in range(ndim)]
        print(u)
        v = derive_by_array(u,x)
        w = derive_by_array(v,x)
        pprint(v)
        pprint(w)
        #(cos(x*t), [x, y, z, t])
        #for key,value in self.dictionary.iteritems():
            #print(key, value)
            #arrays[key.base] = MutableDenseNDimArray.zeros(*key.shape)
            #print(MutableDenseNDimArray([0],key.shape))
        #print "finished"
        #pprint(arrays)
        return
    def get_new_indexedBase(self,shape):
        shape_of_array = tuple([self.ndim for x in range(len(shape))])
        new_Base = IndexedBase('%s%d'%(self.indexedObject_name,self.indexed_object_no), shape= shape_of_array)
        self.indexed_object_no = self.indexed_object_no +1
        print(shape_of_array )
        return new_Base

    def get_functions(self,eq):
        pot = preorder_traversal(eq)
        fns = []
        for p in pot:
            if p in fns:
                pot.skip()
                continue
            elif isinstance(p, EinsteinTerm):
                pot.skip()
                if p.get_indices():
                    fns.append(p)
            elif isinstance(p, Function):
                pot.skip()
                fns.append(p)
            else:
                continue
        return fns
    def get_symbols(self,expression):
        pot = preorder_traversal(expression)
        fns = []
        for p in pot:
            if p in fns:
                pot.skip()
                continue
            elif isinstance(p, EinsteinTerm):
                pot.skip()
                if p.get_indices():
                    fns.append(p.get_indices())
            else:
                continue
        #print(fns)
        return flatten(fns)
    def funtion_to_indexed(self,fn, dictionary):
        print "Non Nested fn is ::  ", fn
        #print(fn)
        indices = []
        if isinstance(fn, EinsteinTerm):
            indices = indices + self.get_symbols(fn)
        elif isinstance(fn, Function):
            for arg in fn.args:
                indices = indices + self.get_symbols(arg)
                #pprint(indices)
        else:
            raise ValueError("function to indexed object failed")
        if indices:
            u = self.get_new_indexedBase(tuple(indices))
            indices = (tuple(indices))
            print "The contraction structure of the function  is  :: ", get_contraction_structure(u[indices])

            print "The indices of the function is  :: ", get_indices(u[indices])
            dictionary[u[indices]] = fn
            print "DICTIONARY ",dictionary
            return u[indices],dictionary
        else:
            return fn,dictionary
        return
    def nested_function_to_indexed(self,fn,dictionary):
        arg_dict = {}
        arg1 = fn.args[0]
        print "Nestedfn is  :: ", fn
        #pprint(fn.args)
        for arg in fn.args[:-1]:
            allfns = self.get_functions(arg)
            #pprint(allfns)
            for ind,fn1 in enumerate(allfns):
                #pprint(fn1)
                temp, arg_dict = self.funtion_to_indexed(fn1,arg_dict)

        #print(arg_dict)
        old_arg = arg1
        for key,value in arg_dict.iteritems():
            if not isinstance(value, EinsteinTerm):
                arg1 = arg1.xreplace({value:key})
        for key,value in arg_dict.iteritems():
            if isinstance(value, EinsteinTerm):
                arg1 = arg1.xreplace({value:key})
        print "The first argument of the nested function is modified to ::  ", arg1
        print "The output structure for the above argument is  :: ", get_indices(arg1)
        print "The contraction structure of the function arg is  :: ", get_contraction_structure(arg1)
        #funtion_to_indexed(fn1,ind,arg_dict)
        indes = list(get_indices(arg1))
        index = list(indes[0])
        u = self.get_new_indexedBase(tuple(index))
        if index:
            arg_dict[u[index]] = arg1
        #index = [u[index]]
        for arg in fn.args[1:]:
            index = index + self.get_symbols(arg)
        u1 = self.get_new_indexedBase(tuple(index))
        out = u1[index[:]]

        nefn = fn.subs(old_arg,arg1)
        print "The input Nested function is updated to :: ",nefn
        arg_dict[out] = nefn
        print "ARG DICT IS", arg_dict
        #dictionary = dictionary.update(arg_dict)
        return out, arg_dict
    def update_dictionary(self,dictionary):
        self.dictionary = dict(self.dictionary.items() + dictionary.items())
        return


class Equation(object):

    """ Describes an equation we want to solve. """

    def __init__(self, expression, ndim, substitutions = [], constants = []):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :returns: None
        """
        local_dict = {'Symbol':EinsteinTerm,'symbols':EinsteinTerm,'Der':Der,'Conservative':Conservative, 'KD':KD}

        self.original = expression

        # Parse the equation.
        self.parsed = parse_expr(self.original, local_dict, evaluate = False)

        # Perform substitutions, if any.
        if substitutions:
            for sub in substitutions:
                temp = parse_expr(sub, local_dict)
                self.parsed = self.parsed.xreplace({temp.lhs: temp.rhs})

        # Update the Einstein Variables in the expression that are constants to is_constant = True
        for term in self.parsed.atoms(EinsteinTerm):
            if any(constant == str(term) for constant in constants):
                term.is_constant = True

        # Expand Einstein terms/indices
        self.expanded = []
        expansion = EinsteinExpansion(self.parsed, ndim)
        LOG.debug("The expanded expression is: %s" % (expansion.expanded))
        # comment out applying formulation
        #if isinstance(expansion.expanded, list):
            #for equation in expansion.expanded:
                #self.expanded.append(self.apply_formulations(equation))
        #else:
            #self.expanded.append(self.apply_formulations(expansion.expanded))

        return

    def apply_formulations(self, expression):
        """ Apply the formulations (e.g. apply Der to create actually differentiate the expanded equations, and KroneckerDelta to actually multiply some terms through by zero). """
        # TODO Move applying Derivative for conservative or Der formulations to their respective classes and call the class from here

        temp_expression = expression
        local_dict = {'Symbol':EinsteinTerm,'symbols':EinsteinTerm,'Der':Der,'Conservative':Conservative}  # TODO: automate from local classes

        if not expression.atoms(Function):
            return expression

        derivative_direction = set()
        # At this point we should not have any other functions except Conservative or Der
        # TODO add a check for the above
        for atom in temp_expression.atoms(Function):
            for arg in atom.args[1:]:
                derivative_direction.add(arg)
        derivative_direction = sorted(derivative_direction, key=lambda x: str(x)[0])  # Sort based on alphabetical order. #FIXME: Currently this only operates on the first letter.
        function_vars = set()
        for term in temp_expression.atoms(EinsteinTerm):
            function_vars.add(term)
        function_vars = function_vars.difference(derivative_direction)
        for term in function_vars:
            if not term.is_constant :
                temp_expression = temp_expression.replace(term, Function('%s'%term)(*derivative_direction))
        new_temp_expression = parse_expr('Eq(%s,%s)'%(str(temp_expression.lhs), str(temp_expression.rhs)),local_dict)
        for atom in new_temp_expression.atoms(Der):
            new_atom = atom.doit()
            new_temp_expression = new_temp_expression.replace(atom,new_atom)
        new_temp_expression = new_temp_expression
        return new_temp_expression


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
