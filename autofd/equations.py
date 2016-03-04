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
from sympy.tensor.index_methods import _get_indices_Mul, _remove_repeated
from sympy.functions.special.tensor_functions import eval_levicivita
import re
import sys
from .array import MutableDenseNDimArray,  derive_by_array, tensorcontraction,tensorproduct,ImmutableDenseNDimArray
import inspect
import sympy.core as core
from sympy import factorial
import numpy as np
import logging
import collections
#from sympy import S, Tuple, MatrixBase
#from sympy import S, Tuple, diff, MatrixBase

#from . import ImmutableDenseNDimArray
from .array import NDimArray
LOG = logging.getLogger(__name__)

# Get Sympy Tensor Functions
#SYMPY_FUNCTIONS = [str(m[0]) for m in inspect.getmembers(tf, inspect.isclass) if m[1].__module__ == tf.__name__]
LOCAL_FUNCTIONS = []


class Conservative(Function):
    LOCAL_FUNCTIONS.append('Conservative')
    @property
    def is_commutative(self):
        return False
    def IndexedObj(self,ndim, indexed_dict, arrays, newarrname):
        # Add repeated calling of the derivatives and functions for support of higher orders
        arguments = {}
        derfn = self.args[0]

        indexobj = IndexedBase('%s'%derfn)
        evaluated, index_struc = evaluate_expression(derfn,arrays,indexed_dict,ndim)
        if not index_struc == None:
            der_struct = index_struc
        else:
            der_struct = []
        if index_struc:
            fns = [indexobj[index_struc]]
        else:
            fns = [indexobj]
        for arg in self.args[1:]:
            der_struct = der_struct+ list(indexed_dict[arg].indices)
            fns.append(indexed_dict[arg])
        shape = []
        for no,ind in enumerate(der_struct):
            if isinstance(ind,Idx):
                shape += [ndim]
            else:
                shape += [ind]
        derivative = MutableDenseNDimArray.zeros(*shape)
        for index in np.ndindex(*shape):
            indexmap = self.split_ind(index,fns)
            derivative[index] = self.apply_der(indexmap,arrays, fns,evaluated)
        outer_indices = remove_repeated_index(der_struct)
        if der_struct:
            derivative = apply_contraction_indexed(outer_indices,der_struct, derivative)
        if outer_indices:
            newouter = [out for out in outer_indices if out!=1]
            if newouter ==outer_indices:
                indobjname = IndexedBase(newarrname,shape=tuple([ndim for x in outer_indices]))[tuple(outer_indices)]
            elif len(newouter) == 1:
                indobjname = IndexedBase(newarrname,shape=tuple([1,ndim]))[tuple(outer_indices)]
                derivative = derivative.reshape(1,ndim)
            indobjname.is_commutative = False
            indexed_dict[self] = indobjname
            arrays[indobjname] = derivative
        else:
            indobjname = EinsteinTerm(newarrname)
            indexed_dict[self] = indobjname
            arrays[indobjname] = derivative
        return arrays, indexed_dict
    def apply_der(self,indexmap,arrays, fns, evaluated):
        if isinstance(fns[0], Indexed):
            derfn = evaluated[indexmap[fns[0]]]
        else:
            derfn = evaluated
        derdire = []
        for fn in fns[1:]:
            if isinstance(fn, Indexed):
                derdire += [arrays[fn][indexmap[fn]]]
            else:
                derdire += [arrays[fn]]

        #derdire = [arrays[fns[i+1]][indexmap[fns[i+1]]] for i,val in enumerate(fns[1:])]
        derivative = Derivative(derfn, *derdire)
        return derivative
    def split_ind(self,index,arrays):
        split = {}
        count = 0
        for arr in arrays:
            if isinstance(arr,Indexed):
                nind = len(arr.indices)
                if nind>0:
                    split[arr] = tuple(index[count:count+nind])
                    count = count + nind
        return split

class KD(Function):
    @property
    def is_commutative(self):
        return False
    # Can combine the two functions below
    def IndexedObj(self,ndim):
        name = str(self.func)
        if len(self.args) >2:
            raise ValueError('Kronecker Delta function should have only two indices')
        ind = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape_of_array = tuple([ndim for x in range(len(ind))])
        indexbase = IndexedBase('%s'%name,shape= shape_of_array)
        indexarray = indexbase[tuple(ind)]
        indexarray.is_commutative=False
        return indexarray
    def ndimarray(self, indexarray):
        array = MutableDenseNDimArray.zeros(*indexarray.shape)
        for index in np.ndindex(*indexarray.shape):
            array[index[:]] = KroneckerDelta(*index)
        return array

class LC(Function):
    LOCAL_FUNCTIONS.append('LC')
    @property
    def is_commutative(self):
        return False
    def IndexedObj(self,ndim):
        name = str(self.func)
        if len(self.args) != 3 or ndim !=3:
            raise ValueError('LeviCivita function should have only three indices indices and third order')
        ind = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape_of_array = tuple([ndim for x in range(len(ind))])
        indexbase = IndexedBase('%s'%name,shape= shape_of_array)
        indexarray = indexbase[tuple(ind)]
        indexarray.is_commutative = False
        return indexarray
    def ndimarray(self, indexarray):
        array = MutableDenseNDimArray.zeros(*indexarray.shape)
        for index in np.ndindex(*indexarray.shape):
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
        return self.indices
    def get_base(self):
        return self.name.split('_')[0]

    def IndexedObj(self,ndim):
        name = self.get_base()
        ind = self.get_indices()
        if len(ind)>1 :
            shape_of_array = tuple([ndim for x in range(len(ind))])
        else:
            ind = [1,self.get_indices()]
            ind = flatten(ind)
            shape_of_array = tuple([1,ndim])
        indexbase = IndexedBase('%s'%name,shape= shape_of_array)
        indexarray = indexbase[tuple(ind)]
        indexarray.is_commutative=False
        return indexarray
    def et_expanded(self, indexmaps):
        newsym = str(self)
        for ind in indexmaps:
            newsym = newsym.replace('_%s'%ind[0], str(ind[1]))
        if self.is_constant:
            new = EinsteinTerm(newsym)
            new.is_constant = True
            return new
        else:
            return EinsteinTerm(newsym)
    def ndimarray(self, indexarray, function=None):
        array = MutableDenseNDimArray.zeros(*indexarray.shape)
        arrayind = indexarray.indices
        for index in np.ndindex(*indexarray.shape):
            indexmap = self.map_indices(arrayind,index)
            val = self.et_expanded(indexmap)
            if not function:
                array[index] = val
            else:
                #array[index] = Function('%s'%val)(*function)
                # CHange here for Function or Indexed Base
                array[index] = IndexedBase('%s'%val)[function]
        return array
    def map_indices(self,arrayind,index):
        maps = []
        for number,ind in enumerate(arrayind):
            if isinstance(ind,Idx):
                maps.append(tuple([ind, index[number]]))
        return maps
def get_index_structure(term):
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
    inds = list(map(get_index_structure, term.args))
    print "MUL indices", inds
    inds = [ind for ind in inds if ind!=None]
    fin_ind = remove_repeated_index(flatten(inds))
    if fin_ind:
        print fin_ind
        print "IN MUL   "
        return fin_ind
    else:
        return None
def get_Add_indices(term):
    """ For additive terms, the indices of the first term is taken as the structure of the additive terms
    """
    inds = list(map(get_index_structure, term.args))
    print(inds)
    if all(ind==None for ind in inds):
        pass
    elif not all([set(x) == set(inds[0]) for x in inds[1:]]):
        raise ValueError("NOT ALL INDICES MATCH in ADD terms of ", term)
    if inds:
        return inds[0]
    else:
        return None
def get_Pow_indices(term):
    base, exp = term.as_base_exp()
    if exp.atoms(Indexed):
        raise ValueError('No indexed objects in exponents  are supported ::', term)
    else:
        base_index_struc = get_index_structure(base)
        if base_index_struc!=None:
            if exp == 2:
                base_index_struc = None
            else:
                raise NotImplementedError("Only Indexed objects to the power 2 are supported")
        else:
            pass
        return base_index_struc
    return
def evaluate_Pow_expression(term, arrays, index_struc):
    base, e = term.as_base_exp()
    print "BASE", base, type(base)
    print "EXPONENT", e, type(e)
    if e.atoms(Indexed):
        raise ValueError('No indexed objects in exponents  are supported ::', term)
    else:
        evaluated, indices = evaluate_Indexed_expression(base, arrays, index_struc)
        print "EVALUATED", evaluated
        if indices!=None:
            if e == 2:
                evaluated = evaluated**e
                tensor_indices = indices
                evaluated = apply_contraction_indexed([], tensor_indices, evaluated)
                indices = []
            else:
                raise NotImplementedError("Only Indexed objects to the power 2 are supported")
        else:
            #print evaluated, "POW", type(evaluated)
            evaluated = Pow((evaluated),(e), evaluate=False)
            #print evaluated, type(evaluated)
            #print srepr(evaluated)
    return evaluated,indices


def evaluate_ADD_expression(term, arrays,index_struc):
    arg_evals = []
    arg_indices = []
    #print "ADD"
    #print "\n"
    for arg in term.args:
        #print "ARG", arg , srepr(arg)
        argeval, arg_index = evaluate_Indexed_expression(arg, arrays, index_struc)
        arg_evals.append(argeval)
        arg_indices.append(arg_index)
        print arg_index
    add_evaluated, indices = add_args(arg_evals,arg_indices)
    return add_evaluated, indices
def add_args(arg_evals,arg_indices):
    if len(arg_evals)==1:
        return arg_evals[0], arg_indices[0]
    # All arguments of addition are scalars, then
    if all([ind==None for ind in arg_indices]):
        evaluated = arg_evals[0]
        for arg in arg_evals[1:]:
            evaluated = evaluated+ arg
        return evaluated, arg_indices[0]
    array_types = (collections.Iterable, MatrixBase, NDimArray)
    print "ADFD"
    print arg_evals
    for number,ind in enumerate(arg_indices):
        if number == 0:
            leading_ind = ind
            evaluated = arg_evals[number]
        else:
            if leading_ind == ind:
                evaluated = evaluated + arg_evals[number]
            else:
                # check the transpose only 2D arrays are supported
                ind.reverse()
                if leading_ind == ind and len(ind) == 2:
                    arr = arg_evals[number]
                    for index in np.ndindex(*arr.shape):
                        transpose = tuple([index[1],index[0]])
                        evaluated[index] = evaluated[index] + arr[transpose]
                else:
                    raise ValueError("IMPLEMNT this in add_args", leading_ind, ind)
    return evaluated, leading_ind

def evaluate_MUL_expression(term, arrays,index_struc):
    out_expression = 1
    tensorprod_indices = []
    for arg in term.args:
        argeval, arg_index = evaluate_Indexed_expression(arg, arrays, index_struc)
        out_expression =tensorproduct(out_expression,argeval)
        if arg_index !=None:
            tensorprod_indices+= arg_index
    indices = remove_repeated_index(tensorprod_indices)
    out_expression = apply_contraction_indexed(indices,tensorprod_indices, out_expression)
    if indices:
        indices = indices
    else:
        indices =None
    return  out_expression, indices

def apply_contraction_indexed(outer_indices, tensor_indices, array):
    contracting_indices = set(tensor_indices).difference(set(outer_indices))
    out_expression = array
    if contracting_indices:
        for index in contracting_indices:
            match = tuple([i for i, x in enumerate(tensor_indices) if x == index])
            out_expression = tensorcontraction(out_expression,match)
            tensor_indices = [i for i  in tensor_indices if i !=index]
    return out_expression
def get_indexed_obj(expression):
    pot = preorder_traversal(expression)
    ETs = []
    Fns = []
    for p in pot:
        if p in ETs+Fns:
            pot.skip()
            continue
        elif isinstance(p, EinsteinTerm):
            pot.skip()
            ETs.append(p)
        elif isinstance(p, Function):
            pot.skip()
            Fns.append(p)
        else:
            continue
    return Fns+ETs

def evaluate_expression(expression, arrays, indexed_dict,ndim):
    indexedobj = get_indexed_obj(expression)
    print '\n'
    print "IN EVALUATE"
    print "INPUT is"
    print expression
    #print arrays
    for Et in indexedobj:
        expression = expression.xreplace({Et:indexed_dict[Et]})
    print "INDEXED IS"
    print expression
    index_struc = get_index_structure(expression)
    print "STRUCTURE is"
    print index_struc
    evaluated, indices = evaluate_Indexed_expression(expression, arrays, index_struc)
    print "Evaluated is"
    print evaluated
    print '\n'
    return evaluated,indices
def evaluate_Indexed_expression(expression, arrays, index_struc):
    # Replace the Einstein terms in the expression by their constituent arrays
    if expression.is_Mul:
        return evaluate_MUL_expression(expression, arrays, index_struc)
    elif expression.is_Add:
        return evaluate_ADD_expression(expression, arrays, index_struc)
        #print evaluated, "ADD"
    elif isinstance(expression, Indexed):
        expindices = list(expression.indices)
        struc = remove_repeated_index(expindices)
        evaluated = apply_contraction_indexed(struc,expindices,arrays[expression])
        return evaluated, struc
    elif isinstance(expression, EinsteinTerm):
        #print arrays[expression], "IN EVALUATE INDEXED", srepr(arrays[expression])
        return arrays[expression],  None
    elif isinstance(expression,Pow):
        return evaluate_Pow_expression(expression, arrays, index_struc)
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
class Der(Function):
    LOCAL_FUNCTIONS.append('Der')
    @property
    def is_commutative(self):
        return False
    def IndexedObj(self,ndim, indexed_dict, arrays, newarrname):
        # Add repeated calling of the derivatives and functions for support of higher orders
        indexobj = IndexedBase('Temp')
        derfn = self.args[0]
        indexobj = IndexedBase('%s'%derfn)
        evaluated, index_struc = evaluate_expression(derfn,arrays,indexed_dict,ndim)
        if not index_struc == None:
            der_struct = index_struc
        else:
            der_struct = []
        if index_struc:
            fns = [indexobj[index_struc]]
        else:
            fns = [indexobj]
        for arg in self.args[1:]:
            der_struct = der_struct+ list(indexed_dict[arg].indices)
            fns.append(indexed_dict[arg])
        shape = []
        for no,ind in enumerate(der_struct):
            if isinstance(ind,Idx):
                shape += [ndim]
            else:
                shape += [ind]
        if shape:
            derivative = MutableDenseNDimArray.zeros(*shape)
            for index in np.ndindex(*shape):
                indexmap = self.split_ind(index,fns)
                derivative[index] = self.apply_der(indexmap,arrays, fns,evaluated)
        else:
            derivative = self.apply_der((0),arrays, fns,evaluated)
        outer_indices = remove_repeated_index(der_struct)
        if der_struct:
            derivative = apply_contraction_indexed(outer_indices,der_struct, derivative)

        if outer_indices:
            newouter = [out for out in outer_indices if out!=1]
            if newouter ==outer_indices:
                indobjname = IndexedBase(newarrname,shape=tuple([ndim for x in outer_indices]))[tuple(outer_indices)]
            elif len(newouter) == 1:
                indobjname = IndexedBase(newarrname,shape=tuple([1,ndim]))[tuple(outer_indices)]
                derivative = derivative.reshape(1,ndim)
            indobjname.is_commutative = False
            indexed_dict[self] = indobjname
            arrays[indobjname] = derivative
        else:
            indobjname = EinsteinTerm(newarrname)
            indexed_dict[self] = indobjname
            arrays[indobjname] = derivative
        return arrays, indexed_dict
    def apply_der(self,indexmap,arrays, fns, evaluated):
        if isinstance(fns[0], Indexed):
            derfn = evaluated[indexmap[fns[0]]]
        else:
            derfn = evaluated
        derdire = []
        for fn in fns[1:]:
            if isinstance(fn, Indexed):
                derdire += [arrays[fn][indexmap[fn]]]
            else:
                derdire += [arrays[fn]]
        derivative = derfn.diff(*derdire)
        return derivative
    def split_ind(self,index,arrays):
        split = {}
        count = 0
        for arr in arrays:
            if isinstance(arr,Indexed):
                nind = len(arr.indices)
                if nind>0:
                    split[arr] = tuple(index[count:count+nind])
                    count = count + nind
        return split

class EinsteinExpansion(object):

    """ Expand an Einstein variable with respect to its indices. """

    def __init__(self, expression, ndim):
        self.expression = expression
        self.expression_indexed = expression
        self.ndim = ndim
        self.expanded = []

        self.is_equality = isinstance(expression, Equality)
        Indexed_dict = {}
        self.indexed_object_no = 0
        self.indexedObject_name = 'Arr'
        ndimarrays = {}
        coord=False
        # update the coordinate ndimarrays:
        for atom in expression.atoms(EinsteinTerm):
            if atom.is_coordinate:
                coord=True
                if atom.get_indices():
                    indict = atom.IndexedObj(self.ndim)
                    arra = atom.ndimarray(indict)
                    coordinates = flatten(arra.tolist())
        # time derivative adding here need to send it out of this
        if coord:
            #print coord
            coordinates = tuple(flatten([coordinates + [EinsteinTerm('t')]]))
            #Indexed_dict[EinsteinTerm('t')] = EinsteinTerm('t')
            #ndimarrays[EinsteinTerm('t')] = EinsteinTerm('t')
        else:
            coordinates = tuple([EinsteinTerm('x0'), EinsteinTerm('x1'),EinsteinTerm('x2'),EinsteinTerm('t')])
        # get the ndim arrays for the Einstein Terms in the equations (u_i,u_j,x_j,x_i) and so on..
        # All the Einstein Terms that are not constants are converted into coordinate functions
        # May be change them later into coordinate indexed objects
        for atom in expression.atoms(EinsteinTerm):
            if atom.is_constant:
                if atom.get_indices():
                    Indexed_dict[atom] = atom.IndexedObj(self.ndim)
                    ndimarrays[Indexed_dict[atom]] = atom.ndimarray(Indexed_dict[atom])
                else:
                    Indexed_dict[atom] = atom
                    ndimarrays[Indexed_dict[atom]] = atom
            else:
                if atom.get_indices() :
                    if atom.get_base() != '':
                        Indexed_dict[atom] = atom.IndexedObj(self.ndim)
                        ndimarrays[Indexed_dict[atom]] = atom.ndimarray(Indexed_dict[atom], coordinates)
                else:
                    Indexed_dict[atom] = atom
                    # CHange here for Function or Indexed Base
                    ndimarrays[Indexed_dict[atom]] = IndexedBase('%s'%atom)[coordinates]

        # get the Ndim arrays for the Kronecker Delta
        for kd in expression.atoms(KD):
            if not kd in Indexed_dict.keys():
                Indexed_dict[kd] = kd.IndexedObj(self.ndim)
                ndimarrays[Indexed_dict[kd]] = kd.ndimarray(Indexed_dict[kd])
        # Do the Functions that are not nested
        for kd in expression.atoms(LC):
            if not kd in Indexed_dict.keys():
                Indexed_dict[kd] = kd.IndexedObj(self.ndim)
                ndimarrays[Indexed_dict[kd]] = kd.ndimarray(Indexed_dict[kd])
        to_eval = []
        for fn in expression.atoms(Function):
            if not fn.args[0].atoms(Function):
                if not fn in Indexed_dict.keys():
                    newarrayname = '%s%d'%(self.indexedObject_name,self.indexed_object_no)
                    self.indexed_object_no = self.indexed_object_no + 1
                    fn.IndexedObj(self.ndim,Indexed_dict,ndimarrays,newarrayname)
            else:
                to_eval.append(fn)
        # evaluate the nested function
        for ev in to_eval:
            #print "in_toeval"
            newarrayname = '%s%d'%(self.indexedObject_name,self.indexed_object_no)
            self.indexed_object_no = self.indexed_object_no + 1
            ev.IndexedObj(self.ndim,Indexed_dict,ndimarrays,newarrayname)
        #now evaluate the RHS of the equation
        print "FINAL"
        evaluated_rhs, rhs_ind = evaluate_expression(expression.rhs, ndimarrays, Indexed_dict,self.ndim)
        evaluated_lhs, lhs_ind = evaluate_expression(expression.lhs, ndimarrays, Indexed_dict,self.ndim)
        newarrayname = 'ALLLHS'
        array_types = (collections.Iterable, MatrixBase, NDimArray)
        finalequation = []
        #pprint(evaluated_rhs)
        if isinstance(evaluated_lhs, array_types):
            print lhs_ind, rhs_ind
            for ind in np.ndindex(evaluated_lhs.shape):
                self.expanded += [Eq(evaluated_lhs[ind], evaluated_rhs[ind])]
                pprint(Eq(evaluated_lhs[ind], evaluated_rhs[ind]))
                pprint(evaluated_lhs[ind])
                pprint(evaluated_rhs[ind])
        else:
            #pprint(evaluated_lhs)
            #pprint(evaluated_rhs)
            self.expanded += [Eq(evaluated_lhs, evaluated_rhs)]
            pprint(Eq(evaluated_lhs, evaluated_rhs))
            #print srepr(evaluated_rhs)

        print "Final expression is", self.expanded
        return


class Equation(object):

    """ Describes an equation we want to solve. """

    def __init__(self, expression, ndim, coordinate_symbol, substitutions = [], constants = []):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :returns: None
        """
        local_dict = {'Symbol':EinsteinTerm,'symbols':EinsteinTerm,'Der':Der,'Conservative':Conservative, 'KD':KD, 'LC':LC}

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
                term.is_commutative = True
            elif term.get_base() == coordinate_symbol or term.get_base() == "t":
                term.is_constant = True
                term.is_coordinate = True


        # Expand Einstein terms/indices

        expansion = EinsteinExpansion(self.parsed, ndim)
        self.expanded = expansion.expanded
        LOG.debug("The expanded expression is: %s" % (expansion.expanded))
        # comment out applying formulation
        #if isinstance(expansion.expanded, list):
            #for equation in expansion.expanded:
                #self.expanded.append(self.apply_formulations(equation))
        #else:
            #self.expanded.append(self.apply_formulations(expansion.expanded))

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
