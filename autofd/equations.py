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
from .array import MutableDenseNDimArray,  derive_by_array, tensorcontraction,tensorproduct
import inspect
import sympy.core as core
from sympy import factorial
import numpy as np
import logging
LOG = logging.getLogger(__name__)

# Get Sympy Tensor Functions
SYMPY_FUNCTIONS = [str(m[0]) for m in inspect.getmembers(tf, inspect.isclass) if m[1].__module__ == tf.__name__]
LOCAL_FUNCTIONS = []


class Conservative(Function):
    LOCAL_FUNCTIONS.append('Conservative')
    @property
    def is_commutative(self):
        return False
    def doit(self,indexed,ndim):
        print "Implement do it for conservative terms"

class KD(Function):
    #LOCAL_FUNCTIONS.append('Der')
    @property
    def is_commutative(self):
        return False
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


def eval_dif(fn, args, **kwargs):
    if kwargs.pop("evaluate", True):
        return diff(fn,*args)
    else:
        return Derivative(fn, *args)

class Der(Function):
    LOCAL_FUNCTIONS.append('Der')
    @property
    def is_commutative(self):
        return False

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

#class IndexedBase(Indexed)

class EinsteinExpansion(object):

    """ Expand an Einstein variable with respect to its indices. """

    def __init__(self, expression, ndim):
        self.expression = expression
        self.expression_indexed = expression
        self.ndim = ndim
        self.expanded = []

        self.is_equality = isinstance(expression, Equality)
        local_functions = self.get_indexed_obj(expression)
        self.dictionary = {}
        self.indexed_object_no = 0
        self.indexedObject_name = 'Arr'
        for ind,fn in enumerate(local_functions):
            dictionary = {}
            if len(fn.atoms(Function)) == 1:
                temp,dictionary = self.funtion_to_indexed(fn,dictionary)
                #self.expression_indexed = self.expression_indexed.xreplace({fn:temp})
                self.update_dictionary(dictionary)
            else:
                temp,dictionary = self.nested_function_to_indexed(fn,dictionary)
                #self.expression_indexed = self.expression_indexed.xreplace({fn:temp})
                self.update_dictionary(dictionary)
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
                    indices = self.index_structure(key)
                    evaluated = self.evaluate_ADD_expression(value,list(indices['outer']), arrays)
                    #addterms[key] = evaluated
                    arrays[key] = evaluated
                    #print arrays[key]
            if muladdterms:
                # reconstruct the term
                for key, value in muladdterms.iteritems():
                    indices = self.index_structure(key)
                    evaluated = self.evaluate_MUL_expression(value,list(indices['outer']), arrays)
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
                index_structure = self.index_structure(term)
                mat = self.evaluate_MUL_expression(term.args,index_struc,arrays)
                add_evaluated = add_evaluated + mat
            else:
                indices = [index for index in term.indices]
                if indices == index_struc:
                    #print "MATCH"
                    #print(arrays[term])
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
                #print arrays[arg], arg
                out_expression = tensorproduct(out_expression,arrays[arg])
                if isinstance(arg,Indexed):
                    tensorprod_indices += [index for index in arg.indices]
                elif isinstance(arg, Expr):
                    #print get_indices(arg)[0]
                    tensorprod_indices += list(get_indices(arg)[0])
                else:
                    raise ValueError('Unknown mul argument')
            else:
                out_expression = tensorproduct(out_expression,arg)
        #print tensorprod_indices
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
        #facs, decomp = self.decompose_args(expression)

        pprint("END")

        return




    def index_structure(self,expression):
        term_list = []
        inner_list = []
        for term, outer, inner in _get_index_structure(expression,
                                                       einstein_notation=True):

            term_list.append(term)
            inner_list.append(inner)

        term_list, inner_list = zip(*sorted(zip(term_list, inner_list),
                                            key=default_sort_key))
        term_list = list(term_list)
        inner_list = list(inner_list)

        return {'monomial_list': term_list, 'outer': outer,
                'inner_list': inner_list}



    def nested_function_to_indexed(self,fn,dictionary):
        arg_dict = {}
        arg1 = fn.args[0]
        old_arg = arg1
        print "Nestedfn is  :: ", fn
        #pprint(fn.args)
        for arg in fn.args[:-1]:
            allfns = self.get_indexed_obj(arg)
            #pprint(allfns)
            for ind,fn1 in enumerate(allfns):
                #pprint(fn1)
                temp, arg_dict = self.funtion_to_indexed(fn1,arg_dict)
                #arg1 = arg1.xreplace({fn1:temp})
        pprint(arg1)
        #print "ARG_DICT", arg_dict
        old_arg = arg1
        #for key,value in arg_dict.iteritems():
            #if not isinstance(value, EinsteinTerm):
                #arg1 = arg1.xreplace({value:key})
        #for key,value in arg_dict.iteritems():
            #if isinstance(value, EinsteinTerm):
                #arg1 = arg1.xreplace({value:key})
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

    def __init__(self, expression, ndim, substitutions = [], constants = []):
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
import collections
def _get_monomial_indices(monomial):
    """List all indices appearing in a monomial with repetition.

    Helper for _get_index_structure."""
    if isinstance(monomial, Mul):
        # Build up list of indices of each term.
        indices_list = []
        for arg in monomial.args:
            indices_list.extend(_get_monomial_indices(arg))
    elif isinstance(monomial, Pow):
        # Get indices of base and duplicate each p times, where p is the
        # power.
        base, power = monomial.args
        indices_list = _get_monomial_indices(base)
        if indices_list and not isinstance(power, Integer):
            msg = "Only integral powers are allowed, not: {!s}"
            raise ValueError(msg.format(power))
        if indices_list:
            indices_list *= power
    elif isinstance(monomial, Indexed):
        indices_list = list(monomial.indices)
    elif isinstance(monomial, Add):
        if monomial.has(Indexed):
            raise ValueError("Not a monomial: {!s}".format(monomial))
        indices_list = []
    else:
        # Assume everything else has no indices.
        indices_list = []
    #indices_list
    #indices_list.sort(key=default_sort_key)
    return indices_list



def _classify_indices(index_list):
    """Classify indices as inner or outer.

    Helper for _get_index_structure.

    Parameters
    ==========

    ``index_list`` : ``list``
        A list of indices appearing in a monomial with repetition.

    Returns
    =======

    tuple : ``set`` of outer indices, ``set`` of inner indices

    """
    outer = set()
    inner = set()
    index_counts = collections.Counter(index_list)
    for index in index_counts:
        if index_counts[index] == 1:
            outer.add(index)
        elif index_counts[index] == 2:
            inner.add(index)
        else:
            msg = "> 2 occurrences of index {!s}: {!s}"
            raise IndexConformanceException(msg.format(index, index_list))
    return outer, inner


def _get_index_structure(expr, einstein_notation):
    """Return a list of tuples describing a tensor expression."""
    expr = expand(expr)
    if isinstance(expr, Add):
        index_structure = []
        for arg in expr.args:
            index_structure.extend(_get_index_structure(arg, einstein_notation))

        # Get outer indices of each term and ensure they're all the same.
        prev_outer = None
        for _, outer, inner in index_structure:
            #print prev_outer, outer, _
            if prev_outer is not None and outer != prev_outer:
                msg = "Inconsistent index structure across sum: {!s}"
                raise IndexConformanceException(msg.format(expr))
            prev_outer = outer

        return index_structure
    else:
        index_list = _get_monomial_indices(expr)
        if einstein_notation:
            outer, inner = _classify_indices(index_list)
        else:
            outer = set(index_list)
            inner = set()
        return [(expr, outer, inner)]
