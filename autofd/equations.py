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
SYMPY_FUNCTIONS = [str(m[0]) for m in inspect.getmembers(tf, inspect.isclass) if m[1].__module__ == tf.__name__]
LOCAL_FUNCTIONS = []


class Conservative(Function):
    LOCAL_FUNCTIONS.append('Conservative')
    @property
    def is_commutative(self):
        return False
    def IndexedObj(self,ndim, indexed_dict, arrays, newarrname):
        # Add repeated calling of the derivatives and functions for support of higher orders
        name = 'Arr_%s%s'
        indexobj = IndexedBase('Temp')
        arguments = {}
        derfn = self.args[0]
        evaluated, index_struc = evaluate_expression(derfn,arrays,indexed_dict,ndim)
        #pprint(arrays)
        arg = self.args[1]
        #print '\n'

        derivative = eval_dif(evaluated,arrays[indexed_dict[arg]], evaluate=False)
        print index_struc,'INDEX STRUCV'

        out_indices = index_struc + list(indexed_dict[arg].indices)
        if out_indices:
            out_struc = list(get_indices(indexobj[tuple(out_indices)])[0])
        all_indices = list(indexed_dict[arg].indices) + index_struc
        #print  all_indices, out_indices
        if len(all_indices)>0:
            der_struc = list(get_indices(indexobj[tuple(all_indices)])[0])
            derivative = apply_contraction_indexed(der_struc,all_indices, derivative)
            #print(derivative)
            if der_struc:
                all_indices =  remove_repeated_index(all_indices)
                out_indices =   remove_repeated_index(out_indices)
                print out_indices, all_indices, der_struc, out_struc
                if out_indices == der_struc:
                    #print "PASS"
                    out_derivative = derivative
                elif len(der_struc) == 2:
                    print "IN HERE"
                    der_struc.reverse()
                    if out_indices == der_struc:
                        mat = transpose(derivative.tomatrix()).as_mutable()
                        mat = MutableDenseNDimArray(mat)
                        out_derivative = mat
                    else:
                        raise ValueError("Index structure is not classified")
                else:
                    raise ValueError("Index structure doesnot match for this")
                indobjname = IndexedBase(newarrname,shape=tuple([ndim for x in out_indices]))[tuple(out_indices)]
                indexed_dict[self] = indobjname
                arrays[indobjname] = out_derivative
            else:
                indobjname = EinsteinTerm(newarrname)
                indexed_dict[self] = indobjname
                arrays[indobjname] = derivative
        else:
            indobjname = EinsteinTerm(newarrname)
            indexed_dict[self] = indobjname
            arrays[indobjname] = derivative
        return arrays, indexed_dict

class KD(Function):
    #LOCAL_FUNCTIONS.append('Der')
    @property
    def is_commutative(self):
        return False
    def IndexedObj(self,ndim):
        name = str(self.func)
        if len(self.args) >2:
            raise ValueError('Kronecker Delta function should have only two indices')
        ind = flatten([p.get_indices() for p in self.args if p.get_indices])
        #print flatten(ind), name
        shape_of_array = tuple([ndim for x in range(len(ind))])
        indexbase = IndexedBase('%s'%name,shape= shape_of_array)
        indexarray = indexbase[tuple(ind)]
        #indexarray.is_commutative=False
        return indexarray
    def ndimarray(self, indexarray):
        array = MutableDenseNDimArray.zeros(*indexarray.shape)
        for index in np.ndindex(*indexarray.shape):
            array[index[:]] = KroneckerDelta(*index)
        return array

    def doit(self,indexed,ndim):

        array = MutableDenseNDimArray.zeros(*indexed.shape)
        if len(indexed.indices) >2:
            raise ValueError('Kronecker Delta function should have only two indices')
        for index in np.ndindex(*indexed.shape):
            array[index[:]] = KroneckerDelta(*index)
        return array


class LeviCivita(Function):
    LOCAL_FUNCTIONS.append('LeviCivita')

    @property
    def is_commutative(self):
        return False

    def doit(self, indexed, ndim):

        n = len(indexed.indices)
        array = MutableDenseNDimArray.zeros(*indexed.shape)
        for index in np.ndindex(*indexed.shape):
            print index
            array[index[:]] = LeviCivita(index)
        #print "in LeviCivita doit",array, self, indexed.indices
        return array


def eval_dif(expr, dx, **kwargs):
    #if kwargs.pop("evaluate", True):
        #return diff(fn,*args)
    #else:
        #return Derivative(fn, *args)

    array_types = (collections.Iterable, MatrixBase, NDimArray)

    if isinstance(dx, array_types):
        dx = MutableDenseNDimArray(dx)
        for i in dx:
            if not i._diff_wrt:
                raise ValueError("cannot derive by this array")

    if isinstance(expr, array_types):
        expr = MutableDenseNDimArray(expr)
        if isinstance(dx, array_types):
            if kwargs.pop("evaluate", True):
                new_array = [[y.diff(x) for y in expr] for x in dx]
            else:
                new_array = [[Derivative(y,x) for y in expr] for x in dx]
            return type(expr)(new_array, dx.shape + expr.shape)
        else:
            if kwargs.pop("evaluate", True):
                return expr.diff(dx)
            else:
                new_array = [Derivative(y,dx) for y in expr]
                return type(expr)(new_array, expr.shape)
    else:
        if isinstance(dx, array_types):
            if kwargs.pop("evaluate", True):
                return MutableDenseNDimArray([expr.diff(i) for i in dx], dx.shape)
            else:
                return MutableDenseNDimArray([Derivative(expr,i) for i in dx], dx.shape)
        else:
            if kwargs.pop("evaluate", True):
                return diff(expr, dx)
            else:
                return Derivative(expr, dx)




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
    def Indexed(self,ndim):
        name = self.name.split('_')[0]
        if name != '':
            shape_of_array= tuple([ndim for x in range(len(self.indices))])
            out = IndexedBase('%s'%name,shape = shape_of_array)
            return out[self.indices]
        else:
            return None
        return
    def get_indices(self):
        return self.indices
    def get_base(self):
        return self.name.split('_')[0]

    def has_index(self, index):
        """ Check to see whether a given index exists in the Einstein variable. """
        return index in self.indices

    def doit(self,indexed,ndim):
        array = MutableDenseNDimArray.zeros(*indexed.shape)
        var = self
        out = []
        # TODO for Einstein term with multiple indices??
        for dim in range(ndim):
            for no,ind in enumerate(self.get_indices()):
                if self.has_index(ind):
                    array[dim]= self.expand_index('_'+str(ind),dim)
        #print "in EISTIEN TERM doit",array, self, indexed.indices
        return array
        #indexed_var = self.Indexed(ndim)
        #out = []
        #indices = self.get_indices()
        #for dim in range(ndim):
            #temp = indexed_var
            #for ind in indices:
                #temp = temp.subs(ind,dim)
            #out.append(temp)
        #print "Expanded EV is ", out
        #return out
    def expand_index(self, index, value):
        """ Replace the index with a particular value (e.g. replace i with 0 in x_i), and return the updated EinsteinTerm object. """
        updated_symbol = str(self).replace(index, str(value))
        if self.is_constant:
            new = EinsteinTerm(updated_symbol)
            new.is_constant = True
            return new
        else:
            return EinsteinTerm(updated_symbol)
    def IndexedObj(self,ndim):
        name = self.get_base()
        ind = self.get_indices()
        if len(ind)>0 :#or self.is_coordinate:
            shape_of_array = tuple([ndim for x in range(len(ind))])
        else:
            ind = [1,self.get_indices()]
            ind = flatten(ind)
            shape_of_array = tuple([1,ndim])
        indexbase = IndexedBase('%s'%name,shape= shape_of_array)
        indexarray = indexbase[tuple(ind)]
        #indexarray.is_commutative=False
        return indexarray
    def et_expanded(self, indexmaps):
        newsym = str(self)
        for ind in indexmaps:
            newsym = newsym.replace('_%s'%ind[0], str(ind[1]))
        #print newsym
        if self.is_constant:
            new = EinsteinTerm(newsym)
            new.is_constant = True
            return new
        else:
            return EinsteinTerm(newsym)
    def ndimarray(self, indexarray, function=None):
        array = MutableDenseNDimArray.zeros(*indexarray.shape)
        #print self, indexarray
        arrayind = indexarray.indices
        for index in np.ndindex(*indexarray.shape):
            indexmap = self.map_indices(arrayind,index)
            val = self.et_expanded(indexmap)
            if not function:
                array[index] = val
            else:
                array[index] = val(*function)
        return array
    def map_indices(self,arrayind,index):
        maps = []
        for number,ind in enumerate(arrayind):
            if isinstance(ind,Idx):
                maps.append(tuple([ind, index[number]]))
        return maps
def evaluate_ADD_expression(terms, index_struc, arrays,ndim):
    if len(index_struc)>0:
        print "in", len(index_struc)
        shape = tuple([ndim for i in index_struc])
        add_evaluated = MutableDenseNDimArray.zeros(*shape)
    else:
        add_evaluated = 0

    for term in terms:
        if term.is_Mul:
            mat, indices = evaluate_MUL_expression(term.args,index_struc,arrays)
            add_evaluated = add_evaluated + mat
        else:
            indices = [index for index in term.indices]
            add_evaluated = add_evaluated + arrays[term]
            # Not required this is handled while evaluating the derivatives and functions
            #if indices == index_struc:
                #add_evaluated = add_evaluated + arrays[term]
            #elif term.rank() == 2:
                #indices.reverse()
                #if indices == index_struc:
                    #pprint(arrays[term].tomatrix())
                    #mat = transpose(arrays[term].tomatrix()).as_mutable()
                    #print "Transpose"
                    #pprint(mat)
                    #mat = MutableDenseNDimArray(mat)
                    #add_evaluated = add_evaluated + mat
            #else:
                #raise ValueError("Indices of %s in add terms dnot match",term)
    #pprint(add_evaluated.tomatrix())
    return add_evaluated
def evaluate_MUL_expression(args,index_struc, arrays):
    #print "IN EVAL MUL"
    out_expression = 1
    tensorprod_indices = []
    for arg in args:
        if arg in arrays.keys():
            out_expression = tensorproduct(out_expression,arrays[arg])
            if isinstance(arg,Indexed):
                tensorprod_indices += [index for index in arg.indices]
            elif isinstance(arg, Expr):
                tensorprod_indices += list(get_indices(arg)[0])
            else:
                raise ValueError('Unknown mul argument')
        else:
            out_expression = tensorproduct(out_expression,arg)
    # get the contracting indices
    print out_expression
    print tensorprod_indices, index_struc
    if tensorprod_indices:
        contracting_indices = set(tensorprod_indices).difference(set(index_struc))
        if contracting_indices:
            for index in contracting_indices:
                match = tuple([i for i, x in enumerate(tensorprod_indices) if x == index])
                out_expression = tensorcontraction(out_expression,match)
                tensorprod_indices = [i for i  in tensorprod_indices if i !=index]
    return  out_expression, tensorprod_indices
#def do_MUL_ADD(expression, ndim_arrays):
    #print "IN MUL ADD"#, expression
    #evaluated = find_index_terms(expression, ndim_arrays)
    ##pprint("END")
    #return evaluated
def find_index_terms(expression,arrays):
    #print "IN FIND TERMS"
    #pprint(expression)
    muladdterms = {}
    addterms = {}
    power_terms = {}
    def _find_terms(expression):
        if expression.is_Mul:
            if expression.atoms(Add):
                muladdterms[expression] = list(expression.args)
                for arg in expression.args:
                    _find_terms(arg)
            else:
                return expression
        elif expression.is_Add:
            terms = []
            for arg in expression.args:
                terms.append(_find_terms(arg))
            addterms[expression] = terms
        elif isinstance(expression, Indexed):
            #indexed.append(expression)
            return expression
            #print "Indexed"
        elif isinstance(expression, EinsteinTerm):
            return expression
        elif expression.is_Atom:
            print "ATOM"
        elif expression.is_Pow or isinstance(expression, exp):
            base, e = expression.as_base_exp()
            if not isinstance(e, Integer):
                raise ValueError('Only integer powers are supported')
            elif not base.atoms(Indexed):
                pass
            else:
                raise NotImplementedError("Implement find term support for Pow")
                #terms = _find_terms(base)
                #powerterms[expression] = terms ** e
        else:
            raise ValueError('I cannot classify the expression',expression)
        return
    _find_terms(expression)
    return muladdterms, addterms,power_terms
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
    # Replace the Einstein terms in the expression by their constituent arrays
    indexedobj = get_indexed_obj(expression)
    print "In evaluate"
    print indexedobj
    pprint(expression)
    for Et in indexedobj:
        expression = expression.xreplace({Et:indexed_dict[Et]})
    pprint(expression)
    if expression.is_Mul and not expression.atoms(Add):
        if not isinstance(expression, Indexed):
            values = list(expression.args)
        else:
            values = [expression]
        #print(_get_indices_Mul(expression, return_dummies=True)); print "IN MUL"
        indices = list(get_indices(expression)[0])
        evaluated, indices = evaluate_MUL_expression(values,indices, arrays)
    elif expression.is_Add:
        print expression, "ADD"
        print expression.atoms(Mul)
        print "EXITINF"
        indices = list(get_indices(expression)[0])
        print indices, expression.args
        evaluated = evaluate_ADD_expression(expression.args,indices, arrays,ndim)
    elif isinstance(expression, Indexed):
        evaluated = arrays[expression]
        indices = list(get_indices(expression)[0])
    elif isinstance(expression, EinsteinTerm):
        evaluated = arrays[expression]
        indices= list(get_indices(expression)[0])
    else:
        muladdterms,addterms,power_terms = find_index_terms(expression,arrays)
        addeval = {}
        if addterms:
            for key, value in addterms.iteritems():
                #if all(arg.atoms(Indexed) for arg in value):
                    indices = get_indices(key)
                    indices = list(get_indices(key)[0])
                    evaluated = evaluate_ADD_expression(value,indices, arrays, ndim)
                    arrays[key] = evaluated
                #elif all(arg.atoms(EinsteinTerm) for arg in value):
                    #indices = []
                    #evaluated = key
                    #for Et in key.atoms(EinsteinTerm):
                        #evaluated =evaluated.subs({Et:arrays[Et]})
                    #arrays[key] = evaluated
                #else:
                    #raise ValueError("Cannot classify ", key,value)
        if muladdterms:
            # reconstruct the term
            for key, value in muladdterms.iteritems():
                indices = list(get_indices(key)[0])
                evaluated,indices = evaluate_MUL_expression(value,indices, arrays)
                arrays[key] = evaluated
    return evaluated,indices
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
        name = 'Arr_%s%s'
        indexobj = IndexedBase('Temp')
        arguments = {}
        derfn = self.args[0]
        evaluated, index_struc = evaluate_expression(derfn,arrays,indexed_dict,ndim)
        #pprint(arrays)
        arg = self.args[1]
        print '\n'

        derivative = eval_dif(evaluated,arrays[indexed_dict[arg]], evaluate=True)
        print index_struc, "INDEX STYRUC"
        index_struc = index_struc
        out_indices = index_struc + list(indexed_dict[arg].indices)
        if out_indices:
            out_struc = list(get_indices(indexobj[tuple(out_indices)])[0])
        all_indices = list(indexed_dict[arg].indices) + index_struc
        #print  all_indices
        if len(all_indices)>0:
            der_struc = list(get_indices(indexobj[tuple(all_indices)])[0])
            derivative = apply_contraction_indexed(der_struc,all_indices, derivative)
            #print(derivative)
            if der_struc:
                all_indices =  remove_repeated_index(all_indices)
                out_indices =   remove_repeated_index(out_indices)
                print self, out_indices, all_indices, der_struc, out_struc
                if out_indices == der_struc:
                    #print "PASS"
                    out_derivative = derivative
                elif len(der_struc) == 2:
                    #print "IN HERE"
                    der_struc.reverse()
                    if out_indices == der_struc:
                        mat = transpose(derivative.tomatrix()).as_mutable()
                        mat = MutableDenseNDimArray(mat)
                        out_derivative = mat
                    else:
                        raise ValueError("Index structure is not classified")
                else:
                    raise ValueError("Index structure doesnot match for this")
                indobjname = IndexedBase(newarrname,shape=tuple([ndim for x in out_indices]))[tuple(out_indices)]
                indexed_dict[self] = indobjname
                arrays[indobjname] = out_derivative
            else:
                indobjname = EinsteinTerm(newarrname)
                indexed_dict[self] = indobjname
                arrays[indobjname] = derivative
        else:
            indobjname = EinsteinTerm(newarrname)
            indexed_dict[self] = indobjname
            arrays[indobjname] = derivative
        return arrays, indexed_dict

    def applyexp(self,indices,values):
        ''' applies the indices to the Einstein terms'''
        out = self
        for no,ind in enumerate(indices):
            for ev in self.atoms(EinsteinTerm):
                if ev.has_index(ind):
                    out = out.replace(ev, ev.expand_index('_'+str(ind),values[no]))
        args = out.args[1:]
        fn = out.args[0]
        return eval_dif(fn,args,evaluate=False)

    def doit(self,indexed,ndim):
        array = MutableDenseNDimArray.zeros(*indexed.shape)
        for index in np.ndindex(*indexed.shape):
            array[index[:]] = self.applyexp(indexed.indices,index)
        #print "in der doit",array, self, indexed.indices
        return array

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
        # update the coordinate ndimarrays:
        #coodinates =
        for atom in expression.atoms(EinsteinTerm):
            if atom.is_coordinate:
                if atom.get_indices():
                    indict = atom.IndexedObj(self.ndim)
                    arra = atom.ndimarray(indict)
                    coordinates = flatten(arra.tolist())
        # time derivative adding here need to send it out of this
        coordinates = tuple(flatten([coordinates + [EinsteinTerm('t')]]))
        ndimarrays[EinsteinTerm('t')] = EinsteinTerm('t')
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
                    ndimarrays[Indexed_dict[atom]] = atom(*coordinates)
        # get the Ndim arrays for the Kronecker Delta and Levicivta functions
        for kd in expression.atoms(KD):
            if not kd in Indexed_dict.keys():
                Indexed_dict[kd] = kd.IndexedObj(self.ndim)
                ndimarrays[Indexed_dict[kd]] = kd.ndimarray(Indexed_dict[kd])
        # do the Levicivita function TODO

        # Do the Functions that are not nested
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
            print "in_toeval"
            newarrayname = '%s%d'%(self.indexedObject_name,self.indexed_object_no)
            self.indexed_object_no = self.indexed_object_no + 1
            ev.IndexedObj(self.ndim,Indexed_dict,ndimarrays,newarrayname)
        #now evaluate the RHS of the equation
        evaluate_expression(expression.rhs, ndimarrays, Indexed_dict,self.ndim)
        #local_functions = self.get_indexed_obj(expression)

        # some useful stuff here
        for key, value in Indexed_dict.iteritems():
            print key, value
            arr= ndimarrays[value]# ndimarrays[value].shape
            pprint(ndimarrays[value])
        #print(ndimarrays)

        #print(to_eval)
        #a1= ndimarrays[Indexed_dict[EinsteinTerm('u_j')]]
        #a2 = ndimarrays[Indexed_dict[EinsteinTerm('x_i')]]
        #a3 = ndimarrays[Indexed_dict[EinsteinTerm('rhou_i')]]
        #io1 = Indexed_dict[EinsteinTerm('u_j')]
        #io2 = Indexed_dict[EinsteinTerm('rhou_i')]
        #tp = tensorproduct(a3,a1)
        #ip = io2*io1
        ##pprint(tp)
        #print '\n\n'
        #print tp.shape, ip, get_indices(ip), get_contraction_structure(ip)
        #tp = tp.reshape(2,2)
        #print tp.shape
        ##pprint(tp.tomatrix())
        #tem = derive_by_array(tp,a2)
        #print tem.shape
        ##tem = tem.reshape(2,2,2)
        ##tem = tensorcontraction(tem,(0,2))
        ##for index in np.ndindex(*tem.shape):
            ##pprint(tem[index])
            ##pprint(index)
        ##print tem.shape
        ##pprint(tem[0])
        #temp2 = derive_by_array(a1,a2)
        #temp2 = temp2.reshape(2,2)
        #print temp2.shape
        #for index in np.ndindex(*temp2.shape):
            #pprint(temp2[index])
            #pprint(index)

        #for ind,fn in enumerate(local_functions):
            ##dictionary = {}
            #if len(fn.atoms(Function)) == 1:
                #print fn
                #temp,dictionary = self.funtion_to_indexed(fn,dictionary)
                ##self.expression_indexed = self.expression_indexed.xreplace({fn:temp})
                #self.update_dictionary(dictionary)
            #else:
                #temp,dictionary = self.nested_function_to_indexed(fn,dictionary)
                ##self.expression_indexed = self.expression_indexed.xreplace({fn:temp})
                #self.update_dictionary(dictionary)
        print "Final expression is", self.expression_indexed

    def get_new_indexedBase(self,shape):
        shape_of_array = tuple([self.ndim for x in range(len(shape))])
        new_Base = IndexedBase('%s%d'%(self.indexedObject_name,self.indexed_object_no), shape= shape_of_array)
        new_Base.is_commutative = False
        self.indexed_object_no = self.indexed_object_no +1
        return new_Base

    def get_indexed_obj(self,eq):
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
        else:
            raise ValueError("function to indexed object failed")
        if indices:
            u = self.get_new_indexedBase(tuple(indices))
            indices = (tuple(indices))
            out = u[indices]
            out.is_commutative= False
            dictionary[fn] = out
            return out,dictionary
        else:
            return fn,dictionary
        return
    def find_index_terms(self,expression,arrays):
        print "IN FIND TERMS"
        pprint(expression)
        muladdterms = {}
        multerms = {}
        inds = []
        addterms = {}
        power_terms = {}
        def _find_terms(expression):
            if expression.is_Mul:
                if expression.atoms(Add):
                    muladdterms[expression] = list(expression.args)
                    for arg in expression.args:
                        _find_terms(arg)
                else:
                    return expression
            elif expression.is_Add:
                terms = []
                for arg in expression.args:
                    terms.append(_find_terms(arg))
                addterms[expression] = terms
            elif isinstance(expression, Indexed):
                #indexed.append(expression)
                return expression
                #print "Indexed"
            elif expression.is_Atom:
                print "ATOM"
            elif expression.is_Pow or isinstance(expression, exp):
                base, e = expression.as_base_exp()
                if not isinstance(e, Integer):
                    raise ValueError('Only integer powers are supported')
                elif not base.atoms(Indexed):
                    pass
                else:
                    raise NotImplementedError("Implement find term support for Pow")
                    #terms = _find_terms(base)
                    #powerterms[expression] = terms ** e
            else:
                raise ValueError('I cannot classify the expression',expression)
            return
        if expression.is_Mul and not expression.atoms(Add):
            terms = [expression]
        elif isinstance(expression, Indexed):
            terms = expression
        else:
            _find_terms(expression)
            pprint(muladdterms)
            pprint(addterms)
            addeval = {}
            if addterms:
                for key, value in addterms.iteritems():
                    indices = get_indices(key)
                    indices = list(get_indices(key)[0])
                    evaluated = self.evaluate_ADD_expression(value,indices, arrays)
                    arrays[key] = evaluated
            if muladdterms:
                # reconstruct the term
                for key, value in muladdterms.iteritems():
                    indices = list(get_indices(key)[0])
                    evaluated = self.evaluate_MUL_expression(value,indices, arrays)
                    arrays[key] = evaluated
                pprint(evaluated)
                print evaluated.shape
                print "MULADD"
        return
    def evaluate_ADD_expression(self,terms, index_struc, arrays):
        print "inADD EVAL"
        shape = tuple([self.ndim for i in index_struc])
        add_evaluated = MutableDenseNDimArray.zeros(*shape)
        for term in terms:
            if term.is_Mul:
                mat = self.evaluate_MUL_expression(term.args,index_struc,arrays)
                add_evaluated = add_evaluated + mat
            else:
                indices = [index for index in term.indices]
                if indices == index_struc:
                    add_evaluated = add_evaluated + arrays[term]
                elif term.rank == 2:
                    indices.reverse()
                    if indices == index_struc:
                        pprint(arrays[term].tomatrix())
                        mat = transpose(arrays[term].tomatrix()).as_mutable()
                        print "Transpose"
                        pprint(mat)
                        mat = MutableDenseNDimArray(mat)
                        add_evaluated = add_evaluated + mat
                else:
                    raise ValueError("Indices of %s in add terms dnot match",term)
        pprint(add_evaluated.tomatrix())
        return add_evaluated
    def evaluate_MUL_expression(self,args,index_struc, arrays):
        print "IN EVAL MUL"
        out_expression = 1
        tensorprod_indices = []
        for arg in args:
            if arg in arrays.keys():
                out_expression = tensorproduct(out_expression,arrays[arg])
                if isinstance(arg,Indexed):
                    tensorprod_indices += [index for index in arg.indices]
                elif isinstance(arg, Expr):
                    tensorprod_indices += list(get_indices(arg)[0])
                else:
                    raise ValueError('Unknown mul argument')
            else:
                out_expression = tensorproduct(out_expression,arg)
        # get the contracting indices
        contracting_indices = set(tensorprod_indices).difference(set(index_struc))
        if contracting_indices:
            for index in contracting_indices:
                match = tuple([i for i, x in enumerate(tensorprod_indices) if x == index])
                out_expression = tensorcontraction(out_expression,match)
                tensorprod_indices = [i for i  in tensorprod_indices if i !=index]
        return  out_expression
    def do_MUL_ADD(self,expression, diction):
        print "IN MUL ADD"#, expression
        for key,value in diction.iteritems():
            if not isinstance(key, EinsteinTerm):
                expression = expression.xreplace({key:value})
        for key,value in diction.iteritems():
            if isinstance(key, EinsteinTerm):
                expression = expression.xreplace({key:value})
        ndim_arrays = {}
        for key1, value1 in diction.iteritems():
            fn = key1
            print fn
            d = get_contraction_structure(value1)
            #pprint(d)
            for key, value in d.iteritems():
                #print value
                if key is not None:
                    #print key, value
                    arr =  fn.func.doit(fn,value1,self.ndim)
                    temp = arr
                    index = value1.indices
                    match = tuple([i for i, x in enumerate(index) if x == key[0]])
                    ndim_arrays[value1] =  temp
                else:
                    temp =  fn.func.doit(fn,value1,self.ndim)
                    ndim_arrays[value1] =  temp
        self.find_index_terms(expression, ndim_arrays)
        pprint("END")
        return

    def nested_function_to_indexed(self,fn,dictionary):
        arg_dict = {}
        arg1 = fn.args[0]
        old_arg = arg1
        print "Nestedfn is  :: ", fn
        for arg in fn.args[:-1]:
            allfns = self.get_indexed_obj(arg)
            for ind,fn1 in enumerate(allfns):
                temp, arg_dict = self.funtion_to_indexed(fn1,arg_dict)
        old_arg = arg1
        indes = list(get_indices(arg1))
        index = list(indes[0])
        u = self.get_new_indexedBase(tuple(index))
        u.is_commutative= False
        self.do_MUL_ADD(arg1,arg_dict)


        if index:
            arg_dict[u[index]] = arg1
        #index = [u[index]]
        for arg in fn.args[1:]:
            index = index + self.get_symbols(arg)
        u1 = self.get_new_indexedBase(tuple(index))
        out = u1[index[:]]
        out.is_commutative= False

        nefn = fn.subs(old_arg,arg1)
        print "The input Nested function is updated to :: ",nefn
        arg_dict[out] = nefn
        #dictionary = dictionary.update(arg_dict)
        return out, arg_dict
    def update_dictionary(self,dictionary):
        self.dictionary = dict(self.dictionary.items() + dictionary.items())
        return


class Equation(object):

    """ Describes an equation we want to solve. """

    def __init__(self, expression, ndim, coordinate_symbol, substitutions = [], constants = []):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :returns: None
        """
        local_dict = {'Symbol':EinsteinTerm,'symbols':EinsteinTerm,'Der':Der,'Conservative':Conservative, 'KD':KD, 'LeviCivita':LeviCivita}

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
            elif term.get_base() == coordinate_symbol or term.get_base() == "t":
                term.is_constant = True
                term.is_coordinate = True


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
        local_dict = {'Symbol':EinsteinTerm,'symbols':EinsteinTerm,'Der':Der,'Conservative':Conservative, 'LeviCivita':LeviCivita}  # TODO: automate from local classes

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
