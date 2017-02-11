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
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, _group_parentheses, _apply_functions, \
    _flatten, _token_splittable
from sympy.tensor.array import MutableDenseNDimArray, tensorcontraction, tensorproduct
from sympy.tensor.array import NDimArray
from sympy.core.function import *
LOG = logging.getLogger(__name__)
from .opensblifunctions import *

LOCAL_FUNCTIONS = []

classdict = {Symbol('Central'): CentralDerivative, Symbol('Temporal'):TemporalDerivative, Symbol('Weno'):WenoDerivative, Symbol('Metric'):MetricDerivative}

class ParsingSchemes(object):
    def set_schemes(cls, expr):
        pot = postorder_traversal(expr)
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
        expr = expr.xreplace(substitutions)
        modifications = [m.xreplace(substitutions) for m in modifications]
        substitutions = {}
        for m in modifications:
            if not (m.args[0].atoms(Derivative).difference(m.args[0].atoms(classdict[cls.args[-1]]))):
                expr = expr.xreplace({m:classdict[cls.args[-1]](m.expr,*m.variables)})
            else:
                raise ValueError("Require better implementation in set schemes")
        return expr

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

class Der(AppliedUndef, ParsingSchemes):

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
        return  expr


class Conservative(AppliedUndef, ParsingSchemes):
    """ Handler for the Conservative function (which uses SymPy's Derivative function). """

    @property
    def is_commutative(cls):
        return False

    def __new__(cls, *args, **kwargs):
        args = set_args(*args, **kwargs)
        ret = super(Conservative, cls).__new__(cls, *args)
        ret.options = kwargs
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
        expr = cls
        local_subs = {}
        for p in pot:
            if isinstance(p, localfuncs):
                #if local_subs:
                    #p = p.subs(local_subs)
                ##print p, p.args, Derivative(p.args[0], p.args[1])
                #pprint([p, p.args[0]])
                #local_subs[p] = Derivative(p.args[0], *p.args[1:-1])
                expr = expr.subs(p, Derivative(p.args[0], *p.args[1:-1]))
                #local_subs += [p]
                #local_subs[p] = Derivative(p.args[0], *p.args[1:-1])
            else:
                continue
        #pprint(local_subs)
        #expr = expr.subs(local_subs)
        # print "Inside conservative"
        # pprint(expr)
        #exit()
        subexpr = expr
        #pprint(expr)
        #pprint(*cls.args[1:-1])
        #pprint(local_subs)
        #pprint(cls)
        #for sub in reversed(local_subs[1:]):
            #expr = expr.subs(sub, evaluate=False)
        #pprint(['jere-1', expr])
        #exit()
        #expr = Derivative(expr, *cls.args[1:-1])
        #pprint(['jere', Derivative(Derivative(cls.args[0], cls.args[1]), cls.args[1])])
        expr = cls.set_schemes(expr)
        return  expr

class MetricDer(AppliedUndef, ParsingSchemes):
    @property
    def is_commutative(cls):
        return False

    def __new__(cls, *args, **kwargs):
        kwargs['scheme'] = 'Central'
        args = set_args(*args, **kwargs)
        ret = super(MetricDer, cls).__new__(cls, *args)
        ret.options = kwargs
        return ret

    def differentiate(cls, Metric):
        localfuncs = (MetricDer)
        notdiff = (ConstantObject, CoordinateObject)
        function_deps  = list(cls.args)[1:-1]
        functions = []
        cls.dummy_index = 0
        mobj = MetricObject('xi_%d'%cls.dummy_index)
        cls.dummy_index = cls.dummy_index +1
        cls.metricobject = mobj
        order = len(cls.atoms(MetricDer))
        if order > 2:
            raise ValueError("Metrics greter than second order derivatives are not tested please uncomment this line and verify metric yourselves")
        pot = postorder_traversal(cls)
        newfn_deps = []
        for p in pot:
            if isinstance(p, localfuncs):
                function_deps += [list(p.args[1:-1])]
                #newfn_deps += [{p:list(p.args[1:-1])}]
            elif isinstance(p, DataObject):
                if p.get_base():
                    functions +=[p]
            else:
                continue
        fndeps = set(flatten(function_deps))
        fns =  set(flatten(functions))
        substits = {}
        revertback = {}
        metric_fns = mobj(*fndeps)
        cls.fndeps = fndeps

        for fn in fns:
            substits[fn] = fn(metric_fns)
        expr = cls.args[0]
        expr = expr.subs(substits)
        pot = postorder_traversal(expr)
        inner_metric_ders = []
        schemes = {}
        for p in pot:
            if isinstance(p, localfuncs):
                inner_metric_ders += [p ]
                newfn_deps += [list(p.args[1:-1])]
                #expr = expr.subs(p, p.args[0].diff(*p.args[1:-1])).doit()
            else:
                continue
        if not inner_metric_ders:
            expr = expr.diff(*cls.args[1:-1])
            expr = expr.doit()
            expr = expr.subs(revertback)
        else:
            if len(inner_metric_ders)> 1:
                raise ValueError("Metrics greter than second order derivatives are not tested please uncomment this line and verify metric yourselves")
            p = inner_metric_ders[0]
            expr = expr.subs(p, p.args[0].diff(*p.args[1:-1])).doit()
            #expr = expr.diff(*cls.args[1:-1])
            expr = expr.doit()
            expr, conversion_subs = cls.convert_metric_ders_dataobjects(expr, order=1)
            expr = expr.diff(*cls.args[1:-1])
            expr = cls.apply_Subs(expr)
            for key, value in conversion_subs.iteritems():
                expr = expr.replace(key,value)
        expr = cls.remove_fn(expr)
        expr = cls.remove_fn(expr)
        expr = cls.set_schemes(expr)
        return expr

    def convert_to_functions(cls, expr):
        """
        Used to convert the expression data objects to functions
        """
        mobj = MetricObject('xi_du%d'%cls.dummy_index)
        cls.dummy_index = cls.dummy_index +1
        metric_fns = mobj(*cls.fndeps)
        for at in expr.atoms(CoordinateObject):
            expr = expr.subs(at, at(metric_fns))
        return expr

    def convert_metric_ders_dataobjects(cls, expression, order):
        # Remove the metric objects functions
        index = 0
        expression = cls.apply_Subs(expression)
        #expression = cls.remove_fn(expression)
        subs = {}
        #expression = cls.apply_Subs(expression)
        for at in expression.atoms(Derivative):
            mobj = MetricObject('xi_%d'%cls.dummy_index)
            cls.dummy_index = cls.dummy_index +1
            metric_fns = mobj(*cls.fndeps)
            f = Function('f%d'%index)
            dobj = CoordinateObject(str(f))
            index = index + 1
            newfn = f(mobj(*cls.fndeps))

            #ppritn(cls.fndeps)
            var_args = []
            for v in at.variables:
                var_args += [cls.remove_fn(v)]
            if isinstance(at.expr, Function):
                if "xi" in str(at.expr.func):
                    #indices = at.expr.get_indices() + [var.get_indices() for var in at.variables]
                    inds = (set(cls.fndeps)).difference(set(var_args))
                    fns = {}
                    for a in at.atoms(Function):
                        arg = var_args + [mobj(*inds)]
                        fns[a] = a.func(*arg)
                    newfn = at.expr.subs(fns)
                    subsfn = Derivative(newfn, *var_args)
                else:
                    fns = {}
                    for a in at.atoms(Function):
                        arg = var_args + [mobj(*cls.fndeps)]
                        fns[a] = a.func(*arg)
                    newfn = at.expr.subs(fns)
                    subsfn = Derivative(newfn, *var_args)
            subs[subsfn] = at
            expression = expression.subs(at, subsfn)
        return expression, subs

    def remove_fn(cls,expr):
        expr = cls.apply_Subs(expr)
        for ind in expr.atoms(Function):
            expr = expr.subs({ind:CoordinateObject(str(ind.func))})
        return expr
    def apply_Subs(cls,expression):
        substis = list(expression.atoms(Subs))
        for work in substis:
            test = list(zip(work.variables, work.point))
            args = work.expr.args
            newargs = [arg.subs(test) for arg in args]
            work1 = type(work.expr)(*newargs)
            expression = expression.subs(work, work1)
        return expression


#from .opensbliequations import *
class Equation(EinsteinStructure):

    """ Describes an equation that is to be solved. """

    def expand(self, expression, ndim, coordinate_symbol, substitutions=[], constants=[], Metric = [], local_dictionary = {}):
        """ Set up an equation, written in Einstein notation, and expand the indices.

        :arg str equation: An equation, written in Einstein notation, and specified in string form.
        :arg int ndim: The dimension of the problem.
        :arg str coordinate_symbol: The spatial coordinate symbol.
        :arg list substitutions: Any substitions to perform (e.g. substituting the stress tensor definition into the Navier-Stokes equations)
        :arg list constants: Any constants like the Reynolds number defined in the equations.
        :returns: None
        """

        #local_dict = {'Symbol': EinsteinTerm, 'symbols': EinsteinTerm, 'Der': Der, 'Conservative': Conservative, 'KD': KD, 'LC': LC, 'Skew': Skew}
        DataSet.dimensions = ndim

        self.original = expression
        constant_dictionary = {}
        for con in constants:
            constant_dictionary[con] = ConstantObject(con)
        local_dict = {'Der': Der, 'Conservative': Conservative, 'KD': KD, 'LC': LC, 'Skew': Skew, 'coordinate':coordinate_symbol,  'time':'t', 'Dot':Dot, 'MetricDer':MetricDer}
        local_dict.update(constant_dictionary)
        #pprint(self.original)
        # Parse the equation.
        self.parsed = parse_expr(self.original, local_dict, tuple([convert_coordinate])+standard_transformations, evaluate=False)

        # Perform substitutions, if any.
        if substitutions:
            for sub in substitutions:
                temp = parse_expr(sub, local_dict, tuple([convert_coordinate])+standard_transformations, evaluate=False)
                self.parsed = self.parsed.subs({temp.lhs: temp.rhs})

        pot = preorder_traversal(self.parsed)
        local = (Der, Conservative)
        for p in pot:
            if isinstance(p, local):
                #pprint([p, p.differentiate()])
                self.parsed = self.parsed.xreplace({p: p.differentiate()})
                pot.skip()
            else:
                continue

        # Handle metric derivatives as special functions
        pot = preorder_traversal(self.parsed)
        local = (MetricDer)
        for p in pot:
            if isinstance(p, local):
                #pprint([p, p.differentiate()])
                self.parsed = self.parsed.xreplace({p: p.differentiate(Metric)})
                pot.skip()
            else:
                continue
        lhs, lhs_indices = self._structure(self.parsed.lhs)
        rhs, rhs_indices = self._structure(self.parsed.rhs)
        lhs = self.substitute_indexed(lhs)
        #pprint(lhs)
        #exit()
        expanded_lhs = self.expand_summations(lhs, ndim)
        rhs = self.substitute_indexed(rhs)
        expanded_rhs = self.expand_summations(rhs, ndim)
        expanded_equation = Eq(expanded_lhs, expanded_rhs)
        if not lhs_indices and not rhs_indices:
            # THIS SHOULD BE MOVED TODO
            expanded_equation = self.apply_functions(expanded_equation)
            #expanded_equation = self.convert_to_data_sets(expanded_equation)
            # Converting to dataObjects
            #for at in expanded_equation.atoms(DataObject):
                #obj = DataSet(str(at)) # By default the location of the dataset is node (0,0,0)
                #expanded_equation = expanded_equation.replace(at, obj)
            ##type_ofeq.add_equations(expanded_equation)
            expanded_equations = expanded_equation

        elif lhs_indices == rhs_indices:
            expanded_equations = expand_free_indices(expanded_equation, lhs_indices, ndim)
            for no, eq in enumerate(expanded_equations):
                eq = self.apply_functions(eq)
                #eq = self.convert_to_data_sets(eq)
                expanded_equations[no] = eq
            #type_ofeq.add_equations(expanded_equations)
        else:
            pprint(lhs_indices)
            pprint(rhs_indices)
            #pprint(expanded_equation)
            raise ValueError("LHS indices and rhs indices are not consistent")
        return expanded_equations


    def apply_functions(self,equation):
        for at in (equation.atoms(Function)):
            if (hasattr(at, 'value')):
                equation= equation.subs(at, at.value)
        return equation
    def convert_to_data_sets(self, equation):
        for at in equation.atoms(DataObject):
            obj = DataSet(str(at)) # By default the location of the dataset is node (0,0,0)
            equation = equation.replace(at, obj)
        return equation



from sympy.parsing.sympy_tokenize import \
    generate_tokens, untokenize, TokenError, \
    NUMBER, STRING, NAME, OP, ENDMARKER

from keyword import iskeyword


def convert_coordinate(tokens, local_dict, global_dict):
    """Trasnforms symbols in parsed expression to OpenSBLI datatypes
    i.e. EinsteinTerm or CoordinateObject or DataObject
    This should take priority over inbuilt sympy transformations
    (standard_transformations) while parsing"""
    result = []
    prevTok = (None, None)

    tokens.append((None, None))  # so zip traverses all tokens
    for tok, nextTok in zip(tokens, tokens[1:]):
        tokNum, tokVal = tok
        nextTokNum, nextTokVal = nextTok
        if tokNum == NAME:
            name = tokVal

            if (local_dict['coordinate']+"_" in name) or local_dict['time'] == name or local_dict['coordinate'] == name:
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

