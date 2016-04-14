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

import logging
import collections

import numpy

from sympy import *
from sympy.parsing.sympy_parser import parse_expr
from sympy.tensor.array import MutableDenseNDimArray, tensorcontraction, tensorproduct
from sympy.tensor.array import NDimArray

LOG = logging.getLogger(__name__)

LOCAL_FUNCTIONS = []


class Skew(Function):

    """ Handler for the energy-conservative formulation of the Navier-Stokes equations, generally referred to as the skew-symmetric formulation.
    This is the Blaisdell version of the skew-symmetric formulation. For more information, see:
    
    [1] G.A. Blaisdell, N.N. Mansour, W.C. Reynolds, Numerical simulations of homogeneous compressible turbulence, Report TF-50, Thermoscience Divison, Department of Mechanical Engineering, Stanford University, Stanford, 1991.
    [2] G.A. Blaisdell, E.T. Spyropoulos, J.H. Qin, The effect of the formulation of nonlinear terms on aliasing errors in spectral methods, Appl. Numer. Math. 1996 207-209

    To get the Blaisdell version of skew symmetric form for the energy and continuity equations, use
        rhou = rho*u    
    and split Conservative((p+rhoE)*u_j,x_j) as Conservative(p*u_j,x_j) + Skew(rhoE*u_j,x_j)
    """

    @property
    def is_commutative(self):
        return False

    def __init__(self, *args):
        return self

    @classmethod
    def eval(cls, *args):
        var = args[0]
        directions = args[1:]
        nvar = len(var.args)
        if nvar == 1:
            return Der(var, *directions)
        elif nvar == 2:
            explicitvars = var.args
            out = 0.5*(Conservative(var, *directions) + explicitvars[1]*Der(explicitvars[0], *directions) + explicitvars[0]*Der(explicitvars[1], *directions))
            return out
        else:
            raise ValueError("More than two terms for Skew formulation.")
        return


class Der(Function):

    """ Handler for the SymPy Derivative function. """

    @property
    def is_commutative(self):
        return False

    def get_indexed(self, ndim, indexed, arrays, new_array_name):

        derivative_function = self.args[0]
        base = IndexedBase('%s' % derivative_function)
        evaluated, index_structure = evaluate_expression(derivative_function, arrays, indexed)
        if index_structure:
            derivative_structure = index_structure
            functions = [base[index_structure]]
        else:
            derivative_structure = []
            functions = [base]

        for arg in self.args[1:]:
            derivative_structure = derivative_structure + list(indexed[arg].indices)
            functions.append(indexed[arg])

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
            for index in numpy.ndindex(*shape):
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
            indexed[self] = indexed_object_name
            arrays[indexed_object_name] = derivative
        else:
            indexed_object_name = EinsteinTerm(new_array_name)
            indexed[self] = indexed_object_name
            arrays[indexed_object_name] = derivative

        return arrays, indexed

    def apply_derivative(self, index_map, arrays, functions, evaluated):
        """ Replace the Derivative calls using SymPy's diff function. """

        if isinstance(functions[0], Indexed):
            function_to_differentiate = evaluated[index_map[functions[0]]]
        else:
            function_to_differentiate = evaluated

        derivative_direction = []
        for fn in functions[1:]:
            if isinstance(fn, Indexed):
                derivative_direction += [arrays[fn][index_map[fn]]]
            else:
                derivative_direction += [arrays[fn]]

        derivative = function_to_differentiate.diff(*derivative_direction)
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

    @property
    def is_commutative(self):
        return False

    def get_indexed(self, ndim, indexed, arrays, new_array_name):

        derivative_function = self.args[0]
        base = IndexedBase('%s' % derivative_function)
        evaluated, index_structure = evaluate_expression(derivative_function, arrays, indexed)

        if index_structure:
            derivative_structure = index_structure
            functions = [base[index_structure]]
        else:
            derivative_structure = []
            functions = [base]

        for arg in self.args[1:]:
            derivative_structure = derivative_structure + list(indexed[arg].indices)
            functions.append(indexed[arg])

        shape = []
        for number, index in enumerate(derivative_structure):
            if isinstance(index, Idx):
                shape += [ndim]
            else:
                shape += [index]

        # Fill an array of the same shape as the derivative structure with the actual SymPy Derivative objects.
        derivative = MutableDenseNDimArray.zeros(*shape)
        for index in numpy.ndindex(*shape):
            index_map = self.split_index(index, functions)
            derivative[index] = self.apply_derivative(index_map, arrays, functions, evaluated)

        # Apply the contraction structure
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
            indexed[self] = indexed_object_name
            arrays[indexed_object_name] = derivative
        else:
            indexed_object_name = EinsteinTerm(new_array_name)
            indexed[self] = indexed_object_name
            arrays[indexed_object_name] = derivative

        return arrays, indexed

    def apply_derivative(self, index_map, arrays, functions, evaluated):
        """ Replace the Conservative calls with Derivative calls (which in turn use SymPy's diff function).

        :returns: The derivative of the function provided, represented as a Derivative object.
        """

        if isinstance(functions[0], Indexed):
            function_to_differentiate = evaluated[index_map[functions[0]]]
        else:
            function_to_differentiate = evaluated

        derivative_direction = []
        for fn in functions[1:]:
            if isinstance(fn, Indexed):
                derivative_direction += [arrays[fn][index_map[fn]]]
            else:
                derivative_direction += [arrays[fn]]

        derivative = Derivative(function_to_differentiate, *derivative_direction)
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

    def get_indexed(self, ndim):
        name = str(self.func)

        if len(self.args) > 2:
            raise ValueError('Kronecker Delta function should have only two indices')

        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape = tuple([ndim for x in range(len(indices))])
        indexed_base = IndexedBase('%s' % name, shape=shape)
        indexed = indexed_base[tuple(indices)]
        indexed.is_commutative = False
        return indexed

    def get_array(self, indexed):
        """ Return an array of KroneckerDelta objects comprising the appropriate indices given in the user's equations. """

        array = MutableDenseNDimArray.zeros(*indexed.shape)
        for index in numpy.ndindex(*indexed.shape):
            array[index[:]] = KroneckerDelta(*index)
        return array


class LC(Function):

    """ Handler for the built-in SymPy LeviCivita function. """

    @property
    def is_commutative(self):
        return False

    def get_indexed(self, ndim):
        name = str(self.func)

        if len(self.args) != 3 or ndim != 3:
            raise ValueError("LeviCivita function should have only three indices.")

        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape = tuple([ndim for x in range(len(indices))])
        indexed_base = IndexedBase('%s' % name, shape=shape)
        indexed = indexed_base[tuple(indices)]
        indexed.is_commutative = False
        return indexed

    def get_array(self, indexed):
        """ Return an array of LeviCivita objects comprising the appropriate indices given in the user's equations. """

        array = MutableDenseNDimArray.zeros(*indexed.shape)
        for index in numpy.ndindex(*indexed.shape):
            array[index[:]] = LeviCivita(*index)
        return array


class EinsteinTerm(Symbol):

    """ Represents any symbol in the equation as a SymPy Symbol object which in turn represents an Einstein term.
    This could be e.g. tau_i_j, but can also be e.g. u_i, rho.
    In other words, all symbols in the equation are Einstein terms, but they can have zero or more indices. """

    is_commutative = False

    def __new__(self, symbol, **assumptions):
        """ Create a new EinsteinTerm.

        :arg str symbol: The symbol under consideration. This can have zero or more indices.
        """

        self._sanitize(assumptions, self)  # Remove any 'None's, etc.
        self.name = str(symbol)

        # Make this into a new SymPy Symbol object.
        self = Symbol.__xnew__(self, self.name, **assumptions)

        # Is this a term that is constant in space and time (i.e. doesn't have any indices).
        self.is_constant = False
        self.is_coordinate = False

        # Extract the indices, which are always preceded by an underscore.
        indices = self.name.split('_')[1:]
        self.indices = [Idx(x) for x in indices]
        return self

    def get_indices(self):
        """ Return a list of the Einstein indices.

        :returns: A list of the Einstein indices.
        :rtype: list
        """
        return self.indices

    def get_base(self):
        """ Return the base name.

        :returns: The base name.
        :rtype: str
        """
        return self.name.split('_')[0]

    def get_indexed(self, ndim):
        """ Given the EinsteinTerm's base name (which gets transformed to an IndexedBase here) and its indices (i.e. Idx objects),
        return an Indexed object.

        :arg int ndim: The dimension of the problem.
        :returns: The EinsteinTerm represented as a SymPy Indexed object.
        """

        name = self.get_base()
        indices = self.get_indices()

        if len(indices) > 0:
            shape = tuple([ndim for x in range(len(indices))])
            indexed_base = IndexedBase('%s' % name, shape=shape)
            indexed = indexed_base[tuple(indices)]
        else:
            indexed = self

        indexed.is_commutative = False
        return indexed

    def get_expanded(self, index_map):
        """ Instantiate a new EinsteinTerm object which is the expanded version of self.

        :arg list index_map: A list of (from, to) tuples which map from each alphabetical/non-expanded Einstein index to the corresponding numerical/expanded one.
        :returns: An expanded version of self, still of the same type (EinsteinTerm)
        :rtype: EinsteinTerm
        """

        expanded = str(self)

        # Expand the indices. This involves replacing the alphabetical Einstein indices with actual numerical indices (e.g. replacing "_j" with "0")
        for index in index_map:
            expanded = expanded.replace('_%s' % index[0], str(index[1]))

        # If the term is a constant then ensure that the relevant flag is set in the new EinsteinTerm object.
        if self.is_constant:
            expanded = EinsteinTerm(expanded)
            expanded.is_constant = True
            # The expanded terms are commutative
            expanded.is_commutative = True
            return expanded
        else:
            expanded = EinsteinTerm(expanded)
            # The expanded terms are commutative
            expanded.is_commutative = True
            return expanded

    def get_array(self, indexed, args=None):
        """ Return an array of Indexed/EinsteinTerm objects.

        :arg sympy.Indexed indexed: A SymPy Indexed term.
        :arg args: A tuple of arguments to be provided to a function (e.g. the LeviCivita function). By default this is None.
        :returns: An array of Indexed/EinsteinTerm objects
        :rtype: sympy.MutableDenseNDimArray, if it has indices else Einstein Term will be returned
        """
        if len(self.get_indices()) == 0:
            value = self
            value.is_commutative = True
            if args:
                array = IndexedBase('%s' % value)[args]
            else:
                array = value
            return array
        array = MutableDenseNDimArray.zeros(*indexed.shape)
        from_indices = indexed.indices
        for index in numpy.ndindex(*indexed.shape):
            index_map = self.map_indices(from_indices, index)
            value = self.get_expanded(index_map)
            value.is_commutative = True
            if args:
                # If this is a function such as the Levi-Civita function, then it will have arguments (e.g. x0, x1, t) which need to be included here.
                array[index] = IndexedBase('%s' % value)[args]
            else:
                array[index] = value

        return array

    def map_indices(self, from_indices, to_indices):
        """ Map each SymPy Idx object to a numerical index.

        :arg tuple from_indices: The tuple of indices to map from.
        :arg tuple to_indices: The tuple of indices to map to.
        :returns: A list of (from, to) tuples of type (Idx, int), e.g. [(j, 1)].
        :rtype: list
        """

        mapping = []
        for number, index in enumerate(from_indices):
            # Only consider indices of type sympy.Idx
            if isinstance(index, Idx):
                mapping.append(tuple([index, to_indices[number]]))
        return mapping


class EinsteinExpansion(object):

    """ Expand an Einstein variable with respect to its indices. """

    def __init__(self, expression, ndim):
        """ Initialise an Einstein expansion system.

        :arg str expression: The expression to expand.
        :arg int ndim: The dimension of the equations.
        """

        self.expression = expression
        self.ndim = ndim
        self.expanded = []

        indexed = {}
        self.indexed_object_number = 0
        self.indexed_object_name = 'Arr'
        arrays = {}
        have_coordinate = False

        # Get the coordinates if the coordinate vector is used in the expression.
        for atom in expression.atoms(EinsteinTerm):
            if atom.is_coordinate:
                have_coordinate = True
                if atom.get_indices():
                    indices = atom.get_indexed(self.ndim)
                    array = atom.get_array(indices)
                    coordinates = flatten(array.tolist())

        # The time symbol is added in here along with the expanded coordinate symbol.
        if have_coordinate:
            coordinates = tuple(flatten([coordinates + [EinsteinTerm('t')]]))
        else:
            # If we do not know the coordinate vector, then generate the coordinates here, and assume that the coordinate vector is called "x".
            coordinates = tuple([EinsteinTerm('x%d' % dim) for dim in range(self.ndim)] + [EinsteinTerm('t')])

        # Get the arrays for the EinsteinTerms in the equations (u_i,u_j,x_j,x_i) and so on.
        # All the Einstein Terms that are not constants are converted into coordinate Indexed Objects.
        for atom in expression.atoms(EinsteinTerm):
            # Constant term
            if atom.is_constant:
                if atom.get_indices():
                    indexed[atom] = atom.get_indexed(self.ndim)
                    arrays[indexed[atom]] = atom.get_array(indexed[atom])
                else:
                    indexed[atom] = atom
                    arrays[indexed[atom]] = atom
            else:
                if atom.get_indices():
                    if atom.get_base():
                        indexed[atom] = atom.get_indexed(self.ndim)
                        arrays[indexed[atom]] = atom.get_array(indexed[atom], coordinates)
                else:
                    indexed[atom] = atom
                    arrays[indexed[atom]] = IndexedBase('%s' % atom)[coordinates]

        # Get the arrays for the Kronecker Delta function
        for kd in expression.atoms(KD):
            if not kd in indexed.keys():
                indexed[kd] = kd.get_indexed(self.ndim)
                arrays[indexed[kd]] = kd.get_array(indexed[kd])

        # Get the arrays for the Levi-Civita function
        for lc in expression.atoms(LC):
            if not lc in indexed.keys():
                indexed[lc] = lc.get_indexed(self.ndim)
                arrays[indexed[lc]] = lc.get_array(indexed[lc])

        # Do the Functions that are not nested
        functions_to_eval = []
        for function in expression.atoms(Function):
            if not function.args[0].atoms(Function):
                if not function in indexed.keys():
                    new_array_name = '%s%d' % (self.indexed_object_name, self.indexed_object_number)
                    self.indexed_object_number = self.indexed_object_number + 1
                    function.get_indexed(self.ndim, indexed, arrays, new_array_name)
            else:
                functions_to_eval.append(function)

        # Evaluate the nested function
        for function in functions_to_eval:
            new_array_name = '%s%d' % (self.indexed_object_name, self.indexed_object_number)
            self.indexed_object_number = self.indexed_object_number + 1
            function.get_indexed(self.ndim, indexed, arrays, new_array_name)

        # Now evaluate the RHS of the equation
        evaluated_rhs, rhs_indices = evaluate_expression(expression.rhs, arrays, indexed)
        evaluated_lhs, lhs_indices = evaluate_expression(expression.lhs, arrays, indexed)

        # Sanity check: Check that the indices on the LHS and RHS match
        if lhs_indices != rhs_indices:
            raise ValueError("Indices of the LHS do not match those of the RHS in the following expression: ", expression)

        # The expanded equations will be of type SymPy Eq.
        array_types = (collections.Iterable, MatrixBase, NDimArray)
        if isinstance(evaluated_lhs, array_types):
            for index in numpy.ndindex(evaluated_lhs.shape):
                self.expanded += [Eq(evaluated_lhs[index], evaluated_rhs[index])]
        else:
            self.expanded += [Eq(evaluated_lhs, evaluated_rhs)]
        return


class Equation(object):

    """ Describes an equation that is to be solved. """

    def __init__(self, expression, ndim, coordinate_symbol, substitutions=[], constants=[]):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :arg int ndim: The dimension of the problem.
        :arg str coordinate_symbol: The spatial coordinate symbol.
        :arg list substitutions: Any substitions to perform (e.g. substituting the stress tensor definition into the Navier-Stokes equations)
        :arg list constants: Any constants like the Reynolds number defined in the equations.
        :returns: None
        """

        local_dict = {'Symbol': EinsteinTerm, 'symbols': EinsteinTerm, 'Der': Der, 'Conservative': Conservative, 'KD': KD, 'LC': LC, 'Skew': Skew}

        self.original = expression

        # Parse the equation.
        self.parsed = parse_expr(self.original, local_dict, evaluate=False)

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
            elif term.get_base() == coordinate_symbol or term.get_base() == "t":  # Spatial and temporal variables are assumed as constant terms here.
                term.is_constant = True
                term.is_coordinate = True
                term.is_commutative = True

        # Expand Einstein terms/indices
        expansion = EinsteinExpansion(self.parsed, ndim)
        self.expanded = expansion.expanded
        # LOG.debug("The expanded expression is: %s" % (expansion.expanded))
        return


def get_index_structure(term):
    """ Get all the Einstein indices of a given term.

    :arg term: The term to get the Einstein indices of.
    :returns: The alphabetical/non-expanded Einstein indices of a given object (or None, if no indices are present).
    """

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


def get_Mul_indices(term):
    """ Get all the Einstein indices in a multiplicative term. """

    # For each of the 'arguments' of a term (e.g. a coefficient like 2/3, or a function like KD[i, j])
    # get the index structure of each one. If no indices exist (e.g. for the coefficient 2/3), then this will be None.
    # For e.g. KD[i, j], this will be [i, j].
    indices = list(map(get_index_structure, term.args))
    # Remove the None's.
    indices = [i for i in indices if i is not None]
    # Remove duplicate indices.
    indices = remove_repeated_index(flatten(indices))
    if indices:
        return indices
    else:
        return None


def get_Add_indices(term):
    """ Get all the Einstein indices in an additive term. The indices of the first term is taken as the structure of the additive terms. """

    indices = list(map(get_index_structure, term.args))
    if all(index is None for index in indices):
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
        if base_index_structure:
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
        raise NotImplementedError("Indexed objects in exponents are not supported: ", term)
    else:
        evaluated, indices = evaluate_Indexed_expression(base, arrays, index_structure)

        if indices:
            if e == 2:
                evaluated = evaluated.applyfunc(lambda i: i**e)
                tensor_indices = indices
                evaluated = apply_contraction([], tensor_indices, evaluated)
                indices = None
            else:
                raise NotImplementedError("Only Indexed objects to the power 2 are supported")
        else:
            evaluated = Pow((evaluated), (e), evaluate=False)
    return evaluated, indices


def evaluate_Add_expression(term, arrays, index_structure):
    """ Evaluate an expression containing an addition operation. """

    arg_evaluations = []
    arg_indices = []
    for arg in term.args:
        arg_evaluated, arg_index = evaluate_Indexed_expression(arg, arrays, index_structure)
        arg_evaluations.append(arg_evaluated)
        arg_indices.append(arg_index)
    add_evaluated, indices = add_args(arg_evaluations, arg_indices)
    return add_evaluated, indices


def add_args(arg_evaluated, arg_indices):
    """ Add all arguments together.

    :arg arg_evaluated: An array of the evaluated EinsteinTerm and/or Function arguments.
    :arg arg_indices: An array of the relevant indices of each of the evaluated EinsteinTerms and/or Functions.
    :returns: The sum of all the evaluated terms.
    """

    # The base case of only one evaluated argument - just return it.
    if len(arg_evaluated) == 1:
        return arg_evaluated[0], arg_indices[0]

    # If all arguments of addition are scalars, then sum them all up and return the result.
    if all([ind is None for ind in arg_indices]):
        evaluated = sum(arg_evaluated)
        return evaluated, arg_indices[0]

    for number, index in enumerate(arg_indices):
        if number == 0:
            # Use the first evaluated argument to initialise the 'evaluated' sum.
            leading_index = index
            evaluated = arg_evaluated[number]
        else:
            if leading_index == index:
                evaluated += arg_evaluated[number]
            else:
                # Check the transpose. Only 2D arrays are supported.
                index.reverse()
                if leading_index == index and len(index) == 2:
                    array = arg_evaluated[number]
                    for index in numpy.ndindex(*array.shape):
                        transpose = tuple([index[1], index[0]])  # Swap the indices around and apply them to the array to get the transpose.
                        evaluated[index] += array[transpose]
                else:
                    raise NotImplementedError("Taking the array transpose is only supported for 2D arrays. ", leading_index, index)

    return evaluated, leading_index


def evaluate_Mul_expression(term, arrays, index_structure):
    """ Evaluate a multiplicative expression. """

    evaluated = 1
    tensor_product_indices = []

    for arg in term.args:
        arg_evaluated, arg_index = evaluate_Indexed_expression(arg, arrays, index_structure)
        evaluated = tensorproduct(evaluated, arg_evaluated)
        if arg_index:
            tensor_product_indices += arg_index

    indices = remove_repeated_index(tensor_product_indices)
    evaluated = apply_contraction(indices, tensor_product_indices, evaluated)
    if indices:
        indices = indices
    else:
        indices = None

    return evaluated, indices


def apply_contraction(outer_indices, tensor_indices, array):
    """ Apply the contraction structure to the array of Indexed objects.

    :arg outer_indices: The outer indices (i.e. any indices that are not repeated/summation indices)
    :arg tensor_indices: The indices that need to be summed over.
    :arg array: The Indexed arrays to apply the index contraction structure to.
    """

    contracting_indices = set(tensor_indices).difference(set(outer_indices))
    result = array
    if contracting_indices:
        for index in contracting_indices:
            match = tuple([i for i, x in enumerate(tensor_indices) if x == index])
            result = tensorcontraction(result, match)
            tensor_indices = [i for i in tensor_indices if i != index]
    return result


def walk_expression_tree(expression):
    """ Perform a pre-order traversal of the expression tree to split up the expression into individual
    EinsteinTerms and functions such as Derivative or Conservative. Any other types (e.g. basic coefficients like 2/3) are skipped.

    :arg expression: The expression from which to build and traverse an expression tree.
    :returns: A list of all individual EinsteinTerms and functions in the expression.
    :rtype: list
    """

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


def evaluate_expression(expression, arrays, indexed):
    """ Evaluate a given expression, thereby expanding the indices and performing any function calls
    (e.g. actually applying the effect of the Kronecker Delta function).

    :arg expression: The expression to evaluate. This will be a SymPy data type.
    :arg dict arrays: A dictionary of (array name, NDimArray) pairs. The NDimArrays contain the expanded variables.
    :arg dict indexed: A dictionary of (non-Indexed, Indexed) pairs. Effectively, each key contains the original
    (non-expanded) EinsteinTerm, and its corresponding value is its equivalent Indexed version. E.g. (rhou_j, rhou[j]).
    :returns: A tuple containing the evaluated expression and a list of indices (e.g. [i,j]).
    :rtype: tuple
    """

    # Get all the EinsteinTerms and functions (e.g. LeviCivita) and replace them by their Indexed versions.
    terms = walk_expression_tree(expression)
    for term in terms:
        expression = expression.xreplace({term: indexed[term]})
    # Find the index structure.
    index_structure = get_index_structure(expression)
    evaluated, indices = evaluate_Indexed_expression(expression, arrays, index_structure)

    return evaluated, indices


def evaluate_Indexed_expression(expression, arrays, index_structure):
    """ Evaluate the Einstein terms in the expression by using their constituent arrays
    and performing the appropriate add/multiply/power operations on them according to the index structure.

    :arg expression: The SymPy expression to evaluate. Each term will be Indexed (as appropriate; constants like "-1" will of course not be Indexed).
    :arg arrays: The arrays containing EinsteinTerms/Indexed components of shape NDim for vectors, or NDim x NDim for rank 2 tensors.
    :arg index_structure: The index structure over which to apply the various mathematical operations.
    """

    if expression.is_Mul:
        return evaluate_Mul_expression(expression, arrays, index_structure)
    elif expression.is_Add:
        return evaluate_Add_expression(expression, arrays, index_structure)
    elif isinstance(expression, Indexed):
        indices = list(expression.indices)
        index_structure = remove_repeated_index(indices)
        evaluated = apply_contraction(index_structure, indices, arrays[expression])
        return evaluated, index_structure
    elif isinstance(expression, EinsteinTerm):
        return arrays[expression], None
    elif isinstance(expression, Pow):
        return evaluate_Pow_expression(expression, arrays, index_structure)
    elif isinstance(expression, Integer):
        return expression, None
    elif expression.is_Atom:
        return expression, None
    else:
        raise ValueError("Unknown expression \"%s\" of type %s" % (expression, type(expression)))
    return


def remove_repeated_index(indices):
    """ Remove duplicate indices in a list of indices.

    :arg list indices: A list of indices to consider.
    :returns: The list of indices with all duplicates removed.
    :rtype: list
    """
    sum_index = {}
    for i in indices:
        if i in sum_index:
            sum_index[i] += 1
        else:
            sum_index[i] = 0
    indices = [x for x in indices if not sum_index[x]]
    return indices


def maximum_derivative_order(equations):
    """ Get the maximum order of the derivatives across the list of equations provided.

    :arg list equations: The list of equations to consider.
    :returns: The maximum order of the derivatives.
    :rtype: int
    """

    order = set()
    for e in equations:
        for atom in e.atoms(Derivative):
            order.add(len(atom.args)-1)
    return max(order)
