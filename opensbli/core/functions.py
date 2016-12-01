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


LOCAL_FUNCTIONS = []

classdict = {Symbol('Central'): CentralDerivative, Symbol('Temporal'):TemporalDerivative, Symbol('Weno'):WenoDerivative}

class EinsteinStructure():
    def structure(cls):
        cls.notation = []
        localfuncs = (Der, KD, Conservative, EinsteinTerm)
        arguments = list(cls.args)[:-1]
        substits = {}
        for no, expr in enumerate(arguments):
            pot = preorder_traversal(expr)
            for p in pot:
                if isinstance(p, localfuncs):
                    if p.structure():
                        substits[p] = p.structure()
                    else:
                        pass
                    pot.skip()
                else:
                    continue
        for no,expr in enumerate(arguments):
            arguments[no] = arguments[no].xreplace(substits)

        cls.der_struct = flatten([get_index_structure(arg) for arg in arguments if get_index_structure(arg)])
        cls.outer_structure = remove_repeated_index(cls.der_struct)
        hashs = str(hash((type(cls).__name__,) + cls.args))
        if cls.outer_structure:
            returning = Indexed("%s%s"%(cls.__class__.__name__, hashs), *cls.outer_structure)
            return  returning
        else:
            return Symbol("%s%s"%(cls.__class__.__name__, hashs))
    def set_schemes(cls, expr):
        pot = postorder_traversal(expr)
        #print("set_schemes", expr.atoms(Derivative))
        #pprint(expr)
        substitutions = {}
        modifications = []
        for d in pot:
            if isinstance(d, Derivative):
                if d.args[0].atoms(Derivative):
                    modifications += [d]
                else:
                    substitutions[d] = classdict[cls.args[-1]](d.expr,*d.variables)
            else:
                continue
        #pprint(srepr(expr))
        #pprint(expr.atoms(CentralDerivative))
        expr = expr.xreplace(substitutions)
        modifications = [m.xreplace(substitutions) for m in modifications]
        substitutions = {}
        for m in modifications:
            if not (m.args[0].atoms(Derivative).difference(m.args[0].atoms(classdict[cls.args[-1]]))):
                expr = expr.xreplace({m:classdict[cls.args[-1]](m.expr,*m.variables)})
            else:
                #pprint(m)
                raise ValueError("Require better implementation in set schemes")
        #pprint(expr)
        #pprint(srepr(expr))

        #for d in expr.atoms(Derivative):
            #scheme = classdict[cls.args[-1]](d.expr,*d.variables)
            #pprint([srepr(d), srepr(classdict[cls.args[-1]](d.expr,*d.variables))])
            #expr = expr.xreplace({d: scheme})

        return expr


class Skew(AppliedUndef):
    """ Handler for the energy-conservative formulation of the Navier-Stokes equations, generally referred to as the skew-symmetric formulation.
    This is the Blaisdell version of the skew-symmetric formulation. For more information, see

    [1] G.A. Blaisdell, N.N. Mansour, W.C. Reynolds, Numerical simulations of homogeneous compressible turbulence, Report TF-50, Thermoscience Divison, Department of Mechanical Engineering, Stanford University, Stanford, 1991.
    [2] G.A. Blaisdell, E.T. Spyropoulos, J.H. Qin, The effect of the formulation of nonlinear terms on aliasing errors in spectral methods, Appl. Numer. Math. 1996 207-209

    To get the Blaisdell version of skew symmetric form for the energy and continuity equations, use rhou = rho*u and split Conservative((p+rhoE)*u_j,x_j) as Conservative(p*u_j,x_j) + Skew(rhoE*u_j,x_j)
    """

    @property
    def is_commutative(cls):
        return False

    def __new__(cls, *args, **kwargs):
        ret = super(Skew, cls).__new__(cls, *args)
        ret.options = kwargs
        pprint(ret)
        ret = ret.applys()
        return ret

    def applys(cls):
        #options = kwargs
        args = cls.args
        kwargs = cls.options
        var = args[0]
        #print var
        directions = args[1:]
        #print directions
        #print kwargs
        nvar = len(var.args)
        if nvar == 1:
            return Der(var, *directions)
        elif nvar == 2:
            explicitvars = var.args
            out = Rational(1,2)*(Conservative(var, *directions, **kwargs) + explicitvars[1]*Der(explicitvars[0], *directions,**kwargs) + explicitvars[0]*Der(explicitvars[1], *directions,**kwargs))
            return out
        else:
            raise ValueError("More than two terms for Skew formulation.")
        return
def set_args(*args, **kwargs):
    args = list(args)
    #pprint([srepr(arg) for arg in args])
    if 'scheme' in kwargs:
        if str(args[-1]) != kwargs['scheme']:
            args = args + [kwargs['scheme']]
    # Include temporal derivatives here
    elif isinstance(args[1], CoordinateObject):
        if args[1].timecoordinate:
            if str(args[-1]) != 'Temporal':
                args = args + ['Temporal']
        else:
            if str(args[-1]) != 'Central':
                args = args + ['Central']
    else:
        if str(args[-1]) != 'Central':
            args = args + ['Central']
    return args


class Der(AppliedUndef, EinsteinStructure):

    """ This is where Der function is applied, as this is a class on its own now
    a. The basic Index structure for each of the derivative and conservative functions
    are processes. Now each Class will have its own index structure evaluation
    b. In simple terms evaluate the Der structure
        i.e Der(u_i,x_j) should be converted to diff(u_i(x_j), x_j)
        and Der(T, x_j,x_j) should be converted to diff(diff(T(x_j),x_j),x_j)
    c. Der(mu*Der(u_i,x_j), x_i) should give
    if mu is constant: then
        mu*Derivative(u_i,x_j,x_i)
    else:
        mu*Derivative(u_i,x_j,x_i) + Derivative(u_i,x_j)* Derivative(mu ,x_i)
    its simple to write an evaluation method

    """

    @property
    def is_commutative(cls):
        return True
    def _print_Der(cls, expr):# NOT working WARNING
        retu = "Derivative(%s)"%(','.join(str(expr.args[:-1])))
        return retu
    def __new__(cls, *args, **kwargs):
        args = set_args(*args, **kwargs)
        ret = super(Der, cls).__new__(cls, *args)
        ret.options = kwargs
        return ret

    def differentiate(cls):
        localfuncs = (Der)
        #pprint(cls)
        #print([a for a in cls.args])
        notdiff = (ConstantObject, CoordinateObject)
        function_deps  = list(cls.args)[1:-1]
        functions = []
        pot = postorder_traversal(cls)
        for p in pot:
            if isinstance(p, localfuncs):
                function_deps += [list(p.args[1:-1])]
            elif isinstance(p, DataObject):
                if p.get_base():
                    functions +=[p]
            else:
                continue
        fndeps = set(flatten(function_deps))
        fns =  set(flatten(functions))
        substits = {}
        revertback = {}
        for fn in fns:
            substits[fn] = fn(*fndeps)
            revertback[fn(*fndeps)] = fn

        #print "PRINTING SUBS"
        #pprint([[type(f) for f in fns], fndeps, substits])
        expr = cls.args[0]
        expr = expr.subs(substits)
        subexpr =expr
        pot = postorder_traversal(expr)
        schemes = {}
        for p in pot:
            if isinstance(p, localfuncs):
                expr = expr.subs(p, p.args[0].diff(*p.args[1:-1])).doit()
            else:
                continue
        expr = expr.diff(*cls.args[1:-1])
        expr = expr.doit()
        expr = expr.subs(revertback)
        expr = cls.set_schemes(expr)

        #pot = preorder_traversal(expr)
        #for d in expr.atoms(Derivative):
            #expr = expr.xreplace({d: CentralDerivative(d.expr,*d.variables)})
        return  expr

    def is_nested(cls):
        if cls.args[0].atoms(Der):
            return False
        else:
            return True

    def to_scheme_ders(cls):
        return



class Conservative(AppliedUndef, EinsteinStructure):
    """ Handler for the Conservative function (which uses SymPy's Derivative function). """

    @property
    def is_commutative(cls):
        return False

    def __new__(cls, *args, **kwargs):
        args = set_args(*args, **kwargs)
        ret = super(Conservative, cls).__new__(cls, *args)
        ret.options = kwargs
        #ret.structure()
        return ret

    def differentiate(cls):
        """
        This requires minimal work as  conservative(rho*u_i,x_j)  ==  Derivative(rho*u_i,x_j)
        TODO Requires more rigorous testing
        """
        localfuncs = (Conservative)
        #notdiff = (ConstantObject, CoordinateObject)
        function_deps  = list(cls.args)[1:-1]
        functions = []
        pot = postorder_traversal(cls)
        expr = cls.args[0]
        for p in pot:
            if isinstance(p, localfuncs):
                expr = expr.subs(p, Derivative(p.args[0], *p.args[1:-1]))
            else:
                continue
        subexpr = expr
        expr = Derivative(expr, *cls.args[1:-1])

        expr = cls.set_schemes(expr)
        #pprint(expr)
        #pprint(srepr(expr))
        #pprint(expr)
        # This the place where we do the CentralDerivative or WenoDerivative
        #pot = preorder_traversal(expr)
        #for d in expr.atoms(Derivative):
            #expr = expr.xreplace({d: CentralDerivative(d.expr,*d.variables)})
        return  expr




class KD(Function):
    """ Handler for the built-in SymPy KroneckerDelta function. """

    @property
    def is_commutative(self):
        return False

    def structure(self):
        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        a = IndexedBase("%s"%self.__class__.__name__)
        a = a[indices]
        a.expression = self
        self.indexed_object = a
        #self.der_struct = a
        return a
    def expansion_symbolic(self):
        return self.indexed_object

    def get_indexed(self, ndim):
        name = str(self.func)

        if len(self.args) > 2:
            raise ValueError('Kronecker Delta function should have only two indices')

        indices = flatten([p.get_indices() for p in self.args if p.get_indices])
        shape = tuple([ndim for x in range(len(indices))])
        indexed_base = DataObject('%s' % name, shape=shape)
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
        indexed_base = DataObject('%s' % name, shape=shape)
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
    In other words, all symbols in the equation are Einstein terms, but they can have zero or more indices.
    Simple way of writing this class would be
    a. Notation, which gives the notaion as an indexed object u[i], tau[i,j] etc. etc.
    b. Value which is nothing but u_i will be u0, u1, u2
    """

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
        #pprint([self, type(self)])

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

    def structure(self):
        if len(self.get_indices()) >0:
            st = IndexedBase("%s"%str(self.get_base()))
            st = st[self.get_indices()]
            st.expression = self
            self.indexed_object = st
            return st
        else:

            st = self
            st.expression = st
            self.indexed_object = st
            #print "here", self
            return st
    def expansion_symbolic(self):
        return self.indexed_object

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
        This is moved to structure
        """
        name = self.get_base()
        indices = self.get_indices()

        if len(indices) > 0:
            shape = tuple([ndim for x in range(len(indices))])
            indexed_base = DataObject('%s' % name, shape=shape)
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
                array = DataObject('%s' % value)[args]
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
                array[index] = DataObject('%s' % value)[args]
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


class ConstantObject(EinsteinTerm):
    is_commutative = True
    def __new__(cls, label):
        ret = super(ConstantObject,cls).__new__(cls, label)
        return ret
    def is_constant(cls):
        return True

class CoordinateObject(EinsteinTerm):
    def __new__(cls, label, **kwargs):
        ret = super(CoordinateObject,cls).__new__(cls, label)
        if 'time' in kwargs:
            ret.timecoordinate = True
        else:
            ret.timecoordinate = False
        return ret

    def get_coordinate_type(cls):
        if cls.timecoordinate:
            return True
        else:
            return False

class DataObject(EinsteinTerm):
    """ This represents the objects this constitutes one of the basic terms if the OpenSBLI objects"""
    def __new__(cls, label, **kw_args):
        ret = super(DataObject, cls).__new__(cls, label, **kw_args)
        return ret

def Sum_to_indexed(expr):
    substits ={}
    for no,s in enumerate(expr.atoms(Sum)):
        sum_indices = set([ind[0] for ind in s.limits])
        indices = (s.function.atoms(Idx).difference(sum_indices))
        if indices:
            temp = Indexed('TEMP%d'%no, *indices)
        else:
            temp = Symbol('TEMP%d'%no)
        substits[temp] = s
        expr = expr.subs(s, temp)
    pprint(expr)
    return expr, substits

class OpenSBLIEquation(Equality):
    def __new__(cls, equation):
        ret = super(OpenSBLIEquation, cls).__new__(cls, equation.lhs, equation.rhs)
        ret.types ={'lhs':0,'rhs':1}
        ret.der_struct = {}
        return ret
    def structure(cls, value, types):
        cls.notation = []
        localfuncs = (KD, DataObject, CoordinateObject, ConstantObject, EinsteinTerm, CentralDerivative, WenoDerivative, TemporalDerivative, LC)
        substits = {}
        # Find the structure of the local functions
        value1 = value
        pot = preorder_traversal(value)
        pprint(value)
        for p in pot:
            if isinstance(p, localfuncs):
                #print(p, type(p))
                if p.structure():
                    substits[p] = p.structure()
                    if isinstance(p, Function):
                        p1 = p.func(*[arg.subs(substits) for arg in p.args])
                    else:
                        p1 = p
                    #p1 = p.subs(substits)
                    #pprint(["in struct eq", p, substits[p], p1])
                    value = value.xreplace({p1: substits[p]})
                    #pprint(value)
                else:
                    pass
                pot.skip()
            elif(isinstance(p, Derivative)):
                pprint(p)
                raise ValueError("This derivative is not classified", p)
            else:
                continue
        #pprint(substits)
        nsubs = {} # Reverse the key and values used for Reverse substitutions
        #for key, v in substits.iteritems():
        


        #pprint(substits)
        #pprint(value.atoms(Pow))
        # Replace the terms in the input expression with their corresponding indexed structures
        for key,val in substits.iteritems():
            #pprint([key, val, val.expression])
            #value = value.subs({key:val})
            nsubs[val] = val.expression
            #pprint(value)
        #value = value.subs(substits, evaluate=False)
        #pprint(value)
        #pprint(value.atoms(Pow))
        # Get the cotaction structure of the expression
        d = (get_contraction_structure(value))
        #pprint(d)
        #from sympy import Expr
        subs = []
        substits_l = {}

        def trav(diction, value):
            outer = []
            replacements = {}
            def mytrav(d1, outer, value, replacements):
                for key in d1:
                    if isinstance(key,Expr):
                        continue
                    for term in d1[key]:
                        if term in d1:
                            for term1 in d1[term]:
                                mytrav(term1, outer, value, replacements)
                    #for v in d1[key]:
                        #v = v.subs
                    if key:
                        for val in d1[key]:
                            summation = val.xreplace(replacements)
                            #summation = val
                            for k in key:
                                summation = Sum(summation, (k,0, Symbol('ndim', integer=True)))
                            #print([val, summation])
                            #value = value.subs({val:summation})
                            if val in replacements:
                                print (exist)
                                raise ValueError("The replacement object already exist this might cause errors")
                            replacements[val] = summation
                        outer += [tuple([key, d1[key]])]

                return replacements, outer
            replacements, outer = mytrav(diction, outer, value, replacements)
            #pprint(outer)
            #pprint(replacements)
            value = value.subs(replacements)
            #pprint(value)
            return value, replacements,outer
        value, replacements, outer = trav(d, value)
        pprint(value)
        def expand_replacements(repl):
            expanded_repl = {}
            for key, value in repl.iteritems():
                #pprint([value, value.variables])
                expanded = 0
                for v in value.variables:
                    for d in range(3):
                        expanded += value.function.xreplace({v:d})
                #pprint([value, expanded])
            return
        expand_replacements(replacements)

        #pprint(nsubs)
        #value= value.xreplace(nsubs)
        #print(value.atoms(Sum))
        #pprint(value.atoms(Sum))
        # Do the traversal of the structure
        #pprint(d.keys())
        #pprint(sorted(d.keys(), key=default_sort_key))
        #for key in d:
            ##pprint(["main ", key])
            #if isinstance(key,Expr):
                #pprint(["EXPR",key, d[key]])
                #for term in d[key]:
                    #pprint(["expr deeper", term])

            #elif key:
                #pprint(["key", key, d[key]])
                ##continue
            ##for term in d[key]:
                ##if term in d:
                    ##pprint(["INTER", key, term])
            ###if key:
                ###pprint(["Final", key,d[key] ])
            ##else:
                ##pprint(["NOIDEA", key, d[key]])
        ##mytrav(d)
        # Find the keys and values in the contraction structure
        #subs, value = _traverse(d, subs, value)
        #
        #pprint(subs)
        #substits_l = {}
        #for s in subs:
            #if s[0]:
                #key, lst = s
                #for v in lst:
                    #v1 = v.xreplace(nsubs)
                    #v1 = v.xreplace(substits_l)
                    #for k in key:
                        #v1 = Sum(v1, (k,0, Symbol('ndim', integer=True)))
                    #substits_l[v] =v1

        #value = value.xreplace(substits_l)
        #pprint(value)
        #value = value.xreplace(nsubs)
        #pprint(value)
        #.xreplace(nsubs)
        return
    def get_lhs_structure(cls):
        return cls.structure(cls.lhs, "lhs")
    def get_rhs_structure(cls):
        return cls.structure(cls.rhs, "rhs")

class Dot(AppliedUndef):

    def __new__(cls, arg1, arg2):
        ret = super(Dot, cls).__new__(cls, arg1, arg2)
        return ret
    def structure(self):
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

        #local_dict = {'Symbol': EinsteinTerm, 'symbols': EinsteinTerm, 'Der': Der, 'Conservative': Conservative, 'KD': KD, 'LC': LC, 'Skew': Skew}

        self.original = expression
        constant_dictionary = {}
        for con in constants:
            constant_dictionary[con] = ConstantObject(con)
        local_dict = {'Der': Der, 'Conservative': Conservative, 'KD': KD, 'LC': LC, 'Skew': Skew, 'coordinate':'x', 'time':'t', 'Dot':Dot}
        local_dict.update(constant_dictionary)
        pprint(self.original)
        # Parse the equation.
        self.parsed = parse_expr(self.original, local_dict, tuple([convert_coordinate])+standard_transformations, evaluate=False)

        # Perform substitutions, if any.
        if substitutions:
            for sub in substitutions:
                temp = parse_expr(sub, local_dict, tuple([convert_coordinate])+standard_transformations, evaluate=False)
                self.parsed = self.parsed.subs({temp.lhs: temp.rhs})
        #pprint(srepr(self.parsed))
        pprint((self.parsed))
        #print(srepr(self.parsed))
        pot = preorder_traversal(self.parsed)
        local = (Der, Conservative)
        for p in pot:
            if isinstance(p, local):
                self.parsed = self.parsed.xreplace({p: p.differentiate()})
                pot.skip()
            else:
                continue
        pprint(srepr(self.parsed))
        pprint(self.parsed)
        # Now obtain the structure of all terms in the parsed equation
        eq = OpenSBLIEquation(self.parsed)
        # Check weather LHS and RHS match
        eq.get_lhs_structure()
        eq.get_rhs_structure()
        #eq.expansion_symbolic()
        print '\n\n'
        return



from sympy.parsing.sympy_tokenize import \
    generate_tokens, untokenize, TokenError, \
    NUMBER, STRING, NAME, OP, ENDMARKER

from keyword import iskeyword


def convert_coordinate(tokens, local_dict, global_dict):
    """Inserts calls to ``Symbol`` for undefined variables."""
    result = []
    prevTok = (None, None)

    tokens.append((None, None))  # so zip traverses all tokens
    for tok, nextTok in zip(tokens, tokens[1:]):
        tokNum, tokVal = tok
        nextTokNum, nextTokVal = nextTok
        if tokNum == NAME:
            name = tokVal

            if (local_dict['coordinate']+"_" in name) or local_dict['time'] == name:
                local_dict['CoordinateObject'] = CoordinateObject
                if local_dict['time'] == name:
                    val = ', **{\'time\':True}'
                    result.extend([
                    (NAME, 'CoordinateObject'),
                    (OP, '('),
                    (NAME, (name)),
                    (OP, ("%s")%val),
                    (OP, ')'),])
                else:
                    result.extend([
                    (NAME, 'CoordinateObject'),
                    (OP, '('),
                    (NAME, (name)),
                    (OP, ')'),])
            elif ("_" in name):
                if name.split("_")[0]:
                    local_dict['DataObject'] = DataObject
                    result.extend([
                        (NAME, 'DataObject'),
                        (OP, '('),
                        (NAME, (name)),
                        (OP, ')'),])
                else:
                    local_dict['EinsteinTerm'] = EinsteinTerm
                    result.extend([
                        (NAME, 'EinsteinTerm'),
                        (OP, '('),
                        (NAME, (name)),
                        (OP, ')'),])
            elif name in local_dict:
                result.append((NAME, name))
            elif name not in global_dict:
                local_dict['DataObject'] = DataObject
                result.extend([
                    (NAME, 'DataObject'),
                    (OP, '('),
                    (NAME, (name)),
                    (OP, ')'),])
            else:
                result.append((NAME, name))
            continue
        else:
            result.append((tokNum, tokVal))
            continue
        prevTok = (tokNum, tokVal)
    return result

