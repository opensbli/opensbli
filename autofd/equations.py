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

import re
import sys
import inspect
import logging
import collections

import numpy as np

from sympy import *
from sympy.core.assumptions import ManagedProperties
from sympy.parsing.sympy_parser import parse_expr
import sympy.functions.special.tensor_functions as tf
from sympy.tensor.index_methods import _get_indices_Mul, _remove_repeated
from sympy.functions.special.tensor_functions import eval_levicivita
import sympy.core as core
from sympy import factorial

from .array import MutableDenseNDimArray, derive_by_array, tensorcontraction,tensorproduct,ImmutableDenseNDimArray
from .array import NDimArray

LOG = logging.getLogger(__name__)

LOCAL_FUNCTIONS = []


class Der(Function):
    
    """ Handler for the SymPy Derivative function. """

    LOCAL_FUNCTIONS.append('Der')
    
    @property
    def is_commutative(self):
        return False
        
    def IndexedObj(self, ndim, indexed_dict, arrays, new_array_name):
        # TODO: Add repeated calling of the derivatives and functions for support of higher orders
        
        derivative_function = self.args[0]
        indexobj = IndexedBase('%s' % derivative_function)
        evaluated, index_structure = evaluate_expression(derivative_function, arrays, indexed_dict)
        if index_structure:
            derivative_structure = index_structure
            functions = [indexobj[index_structure]]
        else:
            derivative_structure = []
            functions = [indexobj]

        for arg in self.args[1:]:
            derivative_structure = derivative_structure + list(indexed_dict[arg].indices)
            functions.append(indexed_dict[arg])

        # Get the shape of the derivative index structure
        shape = []
        for number, index in enumerate(derivative_structure):
            if isinstance(index, Idx):
                shape += [ndim]
            else:
                shape += [index]

        # Apply derivatives (i.e. put SymPy's diff function to use)
        if shape:
            derivative = MutableDenseNDimArray.zeros(*shape)
            for index in np.ndindex(*shape):
                index_map = self.split_index(index, functions)
                derivative[index] = self.apply_derivative(index_map, arrays, functions, evaluated)
        else:
            derivative = self.apply_derivative((0), arrays, functions, evaluated)
            
        # Apply contraction and expand
        outer_indices = remove_repeated_index(derivative_structure)
        if derivative_structure:
            derivative = apply_contraction(outer_indices, derivative_structure, derivative)
        if outer_indices:
            new_outer_indices = [outer for outer in outer_indices if outer != 1]
            if new_outer_indices == outer_indices:
                indexed_object_name = IndexedBase(new_array_name, shape=tuple([ndim for x in outer_indices]))[tuple(outer_indices)]
            else:
                raise ValueError("Indices do not match: ", new_outer_indices, outer_indices)
                
            indexed_object_name.is_commutative = False
            indexed_dict[self] = indexed_object_name
            arrays[indexed_object_name] = derivative
        else:
            indexed_object_name = EinsteinTerm(new_array_name)
            indexed_dict[self] = indexed_object_name
            arrays[indexed_object_name] = derivative
            
        return arrays, indexed_dict
        
    def apply_derivative(self,index_map, arrays, functions, evaluated):
        """ Replace the Derivative calls  """
        if isinstance(functions[0], Indexed):
            derivative_function = evaluated[index_map[functions[0]]]
        else:
            derivative_function = evaluated
        derivative_direction = []
        for fn in functions[1:]:
            if isinstance(fn, Indexed):
                derivative_direction += [arrays[fn][index_map[fn]]]
            else:
                derivative_direction += [arrays[fn]]
        derivative = derivative_function.diff(*derivative_direction)
        return derivative
        
    def split_index(self, index, arrays):
        """ Split up the components of the array. Return a dictionary where each key is an element
        of the array and is associated with the appropriate element of the index tuple. """
        
        split = {}
        count = 0
        for arr in arrays:
            if isinstance(arr, Indexed):
                nind = len(arr.indices)
                if nind > 0:
                    split[arr] = tuple(index[count:count+nind])
                    count = count + nind
        return split


class Conservative(Function):

    """ Handler for the Conservative function (which uses SymPy's Derivative function). """

    LOCAL_FUNCTIONS.append('Conservative')

    @property
    def is_commutative(self):
        return False

    def IndexedObj(self, ndim, indexed_dict, arrays, new_array_name):
        # TODO: Add repeated calling of the derivatives and functions for support of higher orders
        
        arguments = {}
        derivative_function = self.args[0]
        indexobj = IndexedBase('%s' % derivative_function)
        evaluated, index_structure = evaluate_expression(derivative_function, arrays, indexed_dict)
        
        if index_structure:
            derivative_structure = index_structure
            functions = [indexobj[index_structure]]
        else:
            derivative_structure = []
            functions = [indexobj]

        for arg in self.args[1:]:
            derivative_structure = derivative_structure + list(indexed_dict[arg].indices)
            functions.append(indexed_dict[arg])
        
        shape = []
        for number, index in enumerate(derivative_structure):
            if isinstance(index, Idx):
                shape += [ndim]
            else:
                shape += [index]
        
        # Fill an array of the same shape as the derivative structure with the actual SymPy Derivative objects.
        derivative = MutableDenseNDimArray.zeros(*shape)
        for index in np.ndindex(*shape):
            index_map = self.split_index(index,functions)
            derivative[index] = self.apply_derivative(index_map, arrays, functions, evaluated)
        outer_indices = remove_repeated_index(derivative_structure)
        
        # Apply the contraction structure
        if derivative_structure:
            derivative = apply_contraction(outer_indices, derivative_structure, derivative)
        
        if outer_indices:
            new_outer_indices = [outer for outer in outer_indices if outer != 1]
            if new_outer_indices == outer_indices:
                indexed_object_name = IndexedBase(new_array_name, shape=tuple([ndim for x in outer_indices]))[tuple(outer_indices)]
            else:
                raise ValueError("Indices do not match: ", new_outer_indices, outer_indices)
            indexed_object_name.is_commutative = False
            indexed_dict[self] = indexed_object_name
            arrays[indexed_object_name] = derivative
        else:
            indexed_object_name = EinsteinTerm(new_array_name)
            indexed_dict[self] = indexed_object_name
            arrays[indexed_object_name] = derivative
        return arrays, indexed_dict
        
    def apply_derivative(self, index_map, arrays, functions, evaluated):
        if isinstance(functions[0], Indexed):
            derivative_function = evaluated[index_map[functions[0]]]
        else:
            derivative_function = evaluated
        derivative_direction = []
        for fn in functions[1:]:
            if isinstance(fn, Indexed):
                derivative_direction += [arrays[fn][index_map[fn]]]
            else:
                derivative_direction += [arrays[fn]]
        derivative = Derivative(derivative_function, *derivative_direction)
        return derivative
        
    def split_index(self, index, arrays):
        """ Split up the components of the array. Return a dictionary where each key is an element
        of the array and is associated with the appropriate element of the index tuple. """
        
        split = {}
        count = 0
        for arr in arrays:
            if isinstance(arr, Indexed):
                nind = len(arr.indices)
                if nind > 0:
                    split[arr] = tuple(index[count:count+nind])
                    count = count + nind
        return split


class KD(Function):

    """ Handler for the built-in SymPy KroneckerDelta function. """

    @property
    def is_commutative(self):
        return False
        
    # FIXME: Can combine the two functions below
    def IndexedObj(self, ndim):
        name = str(self.func)
        
        if len(self.args) > 2:
            raise ValueError('Kronecker Delta function should have only two indices')
            
        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape_of_array = tuple([ndim for x in range(len(indices))])
        indexed_base = IndexedBase('%s' % name, shape=shape_of_array)
        indexed_array = indexed_base[tuple(indices)]
        indexed_array.is_commutative = False
        return indexed_array
        
    def ndimarray(self, indexed_array):
        """ Return an array of KroneckerDelta objects comprising the appropriate indices given in the user's equations. """
        
        array = MutableDenseNDimArray.zeros(*indexed_array.shape)
        for index in np.ndindex(*indexed_array.shape):
            array[index[:]] = KroneckerDelta(*index)
        return array


class LC(Function):
    
    """ Handler for the built-in SymPy LeviCivita function. """

    LOCAL_FUNCTIONS.append('LC')
    
    @property
    def is_commutative(self):
        return False
        
    def IndexedObj(self, ndim):
        name = str(self.func)
        
        if len(self.args) != 3 or ndim != 3:
            raise ValueError("LeviCivita function should have only three indices.")
            
        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape_of_array = tuple([ndim for x in range(len(indices))])
        indexed_base = IndexedBase('%s' % name, shape=shape_of_array)
        indexed_array = indexed_base[tuple(indices)]
        indexed_array.is_commutative = False
        return indexed_array
        
    def ndimarray(self, indexed_array):
        """ Return an array of LeviCivita objects comprising the appropriate indices given in the user's equations. """
        
        array = MutableDenseNDimArray.zeros(*indexed_array.shape)
        for index in np.ndindex(*indexed_array.shape):
            array[index[:]] = LeviCivita(*index)
        return array


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
        self.is_coordinate = False

        # Extract the indices, which are always preceeded by an underscore.
        indices = self.name.split('_')[1:]
        self.indices = [Idx(x) for x in indices]
        return self
        
    def get_indices(self):
        """ Return a list of the Einstein indices. """
        return self.indices
        
    def get_base(self):
        """ Return the base name. """
        return self.name.split('_')[0]

    def IndexedObj(self,ndim):
        name = self.get_base()
        ind = self.get_indices()
        if len(ind)>0 :
            shape_of_array = tuple([ndim for x in range(len(ind))])
        else:
            ind = [1,self.get_indices()]
            ind = flatten(ind)
            shape_of_array = tuple([1,ndim])
        indexed_base = IndexedBase('%s'%name,shape= shape_of_array)
        indexed_array = indexed_base[tuple(ind)]
        indexed_array.is_commutative=False
        return indexed_array
        
    def get_expanded(self, index_map):
        """ Instantiate a new EinsteinTerm object which is the expanded version of self. """    
    
        expanded = str(self)
        # Expand the indices
        for index in index_map:
            expanded = expanded.replace('_%s' % index[0], str(index[1]))
        # If the term is a constant then ensure that the relevant flag is set in the new EinsteinTerm object.
        if self.is_constant:
            expanded = EinsteinTerm(expanded)
            expanded.is_constant = True
            return expanded
        else:
            return EinsteinTerm(expanded)
            
    def ndimarray(self, indexed_array, function=None):
        """ Return an array of Indexed/EinsteinTerm objects. """
        
        array = MutableDenseNDimArray.zeros(*indexed_array.shape)
        array_indices = indexed_array.indices
        for index in np.ndindex(*indexed_array.shape):
            index_map = self.map_indices(array_indices, index)
            value = self.get_expanded(index_map)
            value.is_commutative = True
            if not function:
                array[index] = value
            else:
                array[index] = IndexedBase('%s' % value)[function]
        return array
        
    def map_indices(self, array_indices, indices):
        """ Map each array index to a SymPy Idx object. """
        
        mapping = []
        for number, index in enumerate(array_indices):
            if isinstance(index, Idx):
                mapping.append(tuple([index, indices[number]]))
        return mapping
       

class EinsteinExpansion(object):

    """ Expand an Einstein variable with respect to its indices. """

    def __init__(self, expression, ndim):
        self.expression = expression
        self.expression_indexed = expression
        self.ndim = ndim
        self.expanded = []

        self.is_equality = isinstance(expression, Equality)
        indexed_dict = {}
        self.indexed_object_number = 0
        self.indexed_object_name = 'Arr'
        ndimarrays = {}
        coord = False
        
        # Update the coordinate ndimarrays
        for atom in expression.atoms(EinsteinTerm):
            if atom.is_coordinate:
                coord = True
                if atom.get_indices():
                    indict = atom.IndexedObj(self.ndim)
                    array = atom.ndimarray(indict)
                    coordinates = flatten(array.tolist())

        # Time derivative added in here along with the expanded coordinate symbol.
        if coord:
            coordinates = tuple(flatten([coordinates + [EinsteinTerm('t')]]))
        else:
            coordinates = tuple([EinsteinTerm('x%d' % dim) for dim in range(self.ndim)] + [EinsteinTerm('t')])

        # Get the ndim arrays for the Einstein Terms in the equations (u_i,u_j,x_j,x_i) and so on.
        # All the Einstein Terms that are not constants are converted into coordinate functions.
        # TODO: Maybe change them later into coordinate indexed objects
        for atom in expression.atoms(EinsteinTerm):
            # Constant term
            if atom.is_constant:
                if atom.get_indices():
                    indexed_dict[atom] = atom.IndexedObj(self.ndim)
                    ndimarrays[indexed_dict[atom]] = atom.ndimarray(indexed_dict[atom])
                else:
                    indexed_dict[atom] = atom
                    ndimarrays[indexed_dict[atom]] = atom
            else:
                if atom.get_indices():
                    if atom.get_base():
                        indexed_dict[atom] = atom.IndexedObj(self.ndim)
                        ndimarrays[indexed_dict[atom]] = atom.ndimarray(indexed_dict[atom], coordinates)
                else:
                    indexed_dict[atom] = atom
                    # Change here for Function or IndexedBase
                    ndimarrays[indexed_dict[atom]] = IndexedBase('%s' % atom)[coordinates]

        # Get the Ndim arrays for the Kronecker Delta function
        for kd in expression.atoms(KD):
            if not kd in indexed_dict.keys():
                indexed_dict[kd] = kd.IndexedObj(self.ndim)
                ndimarrays[indexed_dict[kd]] = kd.ndimarray(indexed_dict[kd])

        # Get the Ndim arrays for the Levi-Civita function
        for lc in expression.atoms(LC):
            if not lc in indexed_dict.keys():
                indexed_dict[lc] = lc.IndexedObj(self.ndim)
                ndimarrays[indexed_dict[lc]] = lc.ndimarray(indexed_dict[lc])

        # Do the Functions that are not nested
        functions_to_eval = []
        for function in expression.atoms(Function):
            if not function.args[0].atoms(Function):
                if not function in indexed_dict.keys():
                    new_array_name = '%s%d' % (self.indexed_object_name, self.indexed_object_number)
                    self.indexed_object_number = self.indexed_object_number + 1
                    function.IndexedObj(self.ndim, indexed_dict, ndimarrays, new_array_name)
            else:
                functions_to_eval.append(function)

        # Evaluate the nested function
        for function in functions_to_eval:
            new_array_name = '%s%d' % (self.indexed_object_name, self.indexed_object_number)
            self.indexed_object_number = self.indexed_object_number + 1
            function.IndexedObj(self.ndim, indexed_dict, ndimarrays, new_array_name)

        # Now evaluate the RHS of the equation
        evaluated_rhs, rhs_ind = evaluate_expression(expression.rhs, ndimarrays, indexed_dict)
        evaluated_lhs, lhs_ind = evaluate_expression(expression.lhs, ndimarrays, indexed_dict)
        
        # Sanity check: Check that the indices on the LHS and RHS match
        if lhs_ind != rhs_ind:
            raise ValueError("Indices of the LHS do not match those of the RHS in the following expression: ", expression)
            
        array_types = (collections.Iterable, MatrixBase, NDimArray)
        if isinstance(evaluated_lhs, array_types):
            for ind in np.ndindex(evaluated_lhs.shape):
                self.expanded += [Eq(evaluated_lhs[ind], evaluated_rhs[ind])]
        else:
            self.expanded += [Eq(evaluated_lhs, evaluated_rhs)]
        return


class Equation(object):

    """ Describes an equation that is to be solved. """

    def __init__(self, expression, ndim, coordinate_symbol, substitutions = [], constants = []):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :arg int ndim: The dimension of the problem.
        :arg str coordinate_symbol: The spatial coordinate symbol.
        :arg list substitutions: Any substitions to perform (e.g. substituting the stress tensor definition into the Navier-Stokes equations)
        :arg list constants: Any constants like the Reynolds number defined in the equations.
        :returns: None
        """

        local_dict = {'Symbol': EinsteinTerm, 'symbols': EinsteinTerm, 'Der': Der, 'Conservative': Conservative, 'KD': KD, 'LC': LC}

        self.original = expression

        # Parse the equation.
        self.parsed = parse_expr(self.original, local_dict, evaluate = False)

        # Perform substitutions, if any.
        if substitutions:
            for sub in substitutions:
                temp = parse_expr(sub, local_dict)
                self.parsed = self.parsed.xreplace({temp.lhs: temp.rhs})

        # Update the Einstein Variables in the expression that are constants
        for term in self.parsed.atoms(EinsteinTerm):
            if any(constant == str(term) for constant in constants):
                term.is_constant = True
                term.is_commutative = True
            elif term.get_base() == coordinate_symbol or term.get_base() == "t": # Spatial and temporal variables are assumed as constant terms here.
                term.is_constant = True
                term.is_coordinate = True
                term.is_commutative = True

        # Expand Einstein terms/indices
        expansion = EinsteinExpansion(self.parsed, ndim)
        self.expanded = expansion.expanded
        LOG.debug("The expanded expression is: %s" % (expansion.expanded))
        
        # TODO: Simplification of the equations

        return


def get_index_structure(term):
    """ Get all the Einstein indices of a given term """

    if isinstance(term, Indexed):
        c = term.indices
        c = remove_repeated_index(list(c))
        return c
    elif term is None:
        return None
    elif term.is_Atom:
        return None
    elif isinstance(term, Idx):
        return [term]
    else:
        if term.is_Mul:
            return get_Mul_indices(term)
        elif term.is_Add:
            return get_Add_indices(term)
        elif term.is_Pow or isinstance(term, exp):
            return get_Pow_indices(term)
    return
    
    
def get_Mul_indices(term):
    """ Get all the Einstein indices in an multiplicative term. """
    
    indices = list(map(get_index_structure, term.args))
    indices = [i for i in indices if i != None]
    # Remove duplicate indices
    indices = remove_repeated_index(flatten(indices))
    if indices:
        return indices
    else:
        return None
        
        
def get_Add_indices(term):
    """ Get all the Einstein indices in an additive term. The indices of the first term is taken as the structure of the additive terms. """
    
    indices = list(map(get_index_structure, term.args))
    if all(ind==None for ind in indices):
        pass
    elif not all([set(x) == set(indices[0]) for x in indices[1:]]):
        raise ValueError("NOT ALL INDICES MATCH in ADD terms of ", term)
    if indices:
        return indices[0]
    else:
        return None
        
        
def get_Pow_indices(term):
    """ Get all the Einstein indices in a Pow term. """

    base, exp = term.as_base_exp()
    if exp.atoms(Indexed):
        raise NotImplementedError('Indexed objects in exponents are not supported: ', term)
    else:
        base_index_structure = get_index_structure(base)
        if base_index_structure!=None:
            if exp == 2:
                base_index_structure = None
            else:
                raise NotImplementedError("Only Indexed objects to the power 2 are supported")
        else:
            pass
        return base_index_structure
    return
    
    
def evaluate_Pow_expression(term, arrays, index_structure):
    """ Evaluate an expression containing an exponent. """

    base, e = term.as_base_exp()
    if e.atoms(Indexed):
        raise NotImplementedError('Indexed objects in exponents are not supported: ', term)
    else:
        evaluated, indices = evaluate_Indexed_expression(base, arrays, index_structure)
        
        if indices:
            if e == 2:
                evaluated = evaluated**e
                tensor_indices = indices
                evaluated = apply_contraction([], tensor_indices, evaluated)
                indices = []
            else:
                raise NotImplementedError("Only Indexed objects to the power 2 are supported")
        else:
            evaluated = Pow((evaluated), (e), evaluate=False)
    return evaluated, indices


def evaluate_Add_expression(term, arrays, index_structure):
    """ Evaluate an expression containing an addition. """

    arg_evals = []
    arg_indices = []
    for arg in term.args:
        arg_eval, arg_index = evaluate_Indexed_expression(arg, arrays, index_structure)
        arg_evals.append(arg_eval)
        arg_indices.append(arg_index)
    add_evaluated, indices = add_args(arg_evals,arg_indices)
    return add_evaluated, indices
    
    
def add_args(arg_evals, arg_indices):
    """ Add all arguments together. """

    if len(arg_evals)==1:
        return arg_evals[0], arg_indices[0]
        
    # If all arguments of addition are scalars, then:
    if all([ind==None for ind in arg_indices]):
        evaluated = arg_evals[0]
        for arg in arg_evals[1:]:
            evaluated = evaluated+ arg
        return evaluated, arg_indices[0]
        
    array_types = (collections.Iterable, MatrixBase, NDimArray)
    
    for number, ind in enumerate(arg_indices):
        if number == 0:
            leading_index = ind
            evaluated = arg_evals[number]
        else:
            if leading_index == ind:
                evaluated = evaluated + arg_evals[number]
            else:
                # Check the transpose. Only 2D arrays are supported.
                ind.reverse()
                if leading_index == ind and len(ind) == 2:
                    arr = arg_evals[number]
                    for index in np.ndindex(*arr.shape):
                        transpose = tuple([index[1],index[0]])
                        evaluated[index] = evaluated[index] + arr[transpose]
                else:
                    raise NotImplementedError("Taking the array transpose is only supported for 2D arrays.", leading_index, ind)
                    
    return evaluated, leading_index


def evaluate_Mul_expression(term, arrays, index_structure):
    """ Evaluate a multiplicative expression. """

    evaluated = 1
    tensorprod_indices = []
    
    for arg in term.args:
        arg_eval, arg_index = evaluate_Indexed_expression(arg, arrays, index_structure)
        evaluated = tensorproduct(evaluated, arg_eval)
        if arg_index:
            tensorprod_indices += arg_index
            
    indices = remove_repeated_index(tensorprod_indices)
    evaluated = apply_contraction(indices, tensorprod_indices, evaluated)
    if indices:
        indices = indices
    else:
        indices = None
        
    return evaluated, indices


def apply_contraction(outer_indices, tensor_indices, array):
    """ Apply the contraction structure to the array of Indexed objects. """

    contracting_indices = set(tensor_indices).difference(set(outer_indices))
    result = array
    if contracting_indices:
        for index in contracting_indices:
            match = tuple([i for i, x in enumerate(tensor_indices) if x == index])
            result = tensorcontraction(result, match)
            tensor_indices = [i for i in tensor_indices if i != index]
            
    return result
    
    
def get_indexed_object(expression):
    """ Perform a pre-order traversal of the expression tree to get all Indexed
    EinsteinTerms and functions such as Derivative or Conservative. """
    
    pot = preorder_traversal(expression)
    einstein_terms = []
    functions = []
    for p in pot:
        if p in einstein_terms+functions:
            pot.skip()
            continue
        elif isinstance(p, EinsteinTerm):
            pot.skip()
            einstein_terms.append(p)
        elif isinstance(p, Function):
            pot.skip()
            functions.append(p)
        else:
            continue
    return functions + einstein_terms


def evaluate_expression(expression, arrays, indexed_dict):
    indexed_object = get_indexed_object(expression)
    for einstein_term in indexed_object:
        expression = expression.xreplace({einstein_term:indexed_dict[einstein_term]})
    index_structure = get_index_structure(expression)
    evaluated, indices = evaluate_Indexed_expression(expression, arrays, index_structure)
    return evaluated, indices
    
    
def evaluate_Indexed_expression(expression, arrays, index_structure):
    """ Replace the Einstein terms in the expression by their constituent arrays. """
    
    if expression.is_Mul:
        return evaluate_Mul_expression(expression, arrays, index_structure)
    elif expression.is_Add:
        return evaluate_Add_expression(expression, arrays, index_structure)
    elif isinstance(expression, Indexed):
        expindices = list(expression.indices)
        struc = remove_repeated_index(expindices)
        evaluated = apply_contraction(struc,expindices,arrays[expression])
        return evaluated, struc
    elif isinstance(expression, EinsteinTerm):
        return arrays[expression],  None
    elif isinstance(expression,Pow):
        return evaluate_Pow_expression(expression, arrays, index_structure)
    elif isinstance(expression, Integer):
        return expression, None
    elif expression.is_Atom:
        return expression, None
    else:
        raise ValueError("SOME THING NOT UNDERSTOOD", expression, type(expression))
    return
    
    
def remove_repeated_index(listofind):
    sum_index = {}
    for i in listofind:
        if i in sum_index:
            sum_index[i] += 1
        else:
            sum_index[i] = 0
    listofind = [x for x in listofind if not sum_index[x]]
    return listofind
