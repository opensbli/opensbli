#!/usr/bin/env python
#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, David J. Lusher, Neil D. Sandham. 

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
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>

import sys
from math import ceil

# Import local utility functions
from sympy.tensor.array import MutableDenseNDimArray, tensorcontraction
from opensbli.core.opensbliobjects import * 
from sympy import * 
from opensbli.core.opensbliobjects import ConstantObject
from opensbli.core.opensblifunctions import CentralDerivative
from opensbli.core.opensbliequations import OpenSBLIExpression
from sympy.parsing.sympy_parser import parse_expr

class EulerEquations(object):
    def __init__(self, ndim, **kwargs):
        """ """
        # self.constants = []
        # self.general_equations()
        # self.formulas = []
        # self.velocity_components()
        # self.pressure_equation()
        # self.speed_of_sound()
        self.ndim = ndim
        # coordinate_symbol = "x"
        # metrics = []
        # problem = Problem(self.equations, [], self.ndim, self.constants, coordinate_symbol, self.formulas)
        # self.expanded_equations = problem.get_expanded(problem.equations)
        # self.expanded_formulas = problem.get_expanded(problem.formulas)
        # self.problem = problem
        # self.time_derivative = EinsteinTerm('t')
        # self.space_derivative_symbols = EinsteinTerm('%s_j'%(coordinate_symbol))
        # self.space_derivative_symbols= self.space_derivative_symbols.get_array(self.space_derivative_symbols.get_indexed(self.ndim))
        # self.eq_to_vector_form()
        # self.formula_to_dict()
        self.REV = {}
        self.LEV = {}
        self.EV= {}
        # self.conser_to_primitive()
        return

    def generate_eig_system(self, block=None):
        ndim = self.ndim
        # coordinates = block.coordinates
        coordinate_symbol = "x"
        cart = CoordinateObject('%s_i'%(coordinate_symbol))
        coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        pprint(coordinates)
        # # Check if block has metrics to be used, else revert to cartesian
        # if block.metric_transformations:
        #     met_symbols = block.FD_simple_symbols
        #     pprint(block.FD_eval_symbols
        met_symbols = eye(ndim)

        if ndim == 1:
            list1 = ['H']
            list2 = ['a**2/(gama-1) + u0**2/2']
            list1 = [parse_expr(l) for l in list1]
            list2 = [parse_expr(l) for l in list2]

            ev = 'diag([u0-a, u0, u0+a])'
            REV = 'Matrix([[1,1,1], [u0-a,u0,u0+a], [H-u0*a,u0**2 /2,H+u0*a]])'
            LEV = 'Matrix([[ u0*(2*H + a*u0 - u0**2)/(2*a*(2*H - u0**2)), (-H - a*u0 + u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)],[2*(H - u0**2)/(2*H - u0**2),2*u0/(2*H - u0**2), 2/(-2*H + u0**2)],[u0*(-2*H + a*u0 + u0**2)/(2*a*(2*H - u0**2)),  (H - a*u0 - u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)]])'
            ev = parse_expr(ev)
            REV = parse_expr(REV)
            LEV = parse_expr(LEV)

            subs_dict = dict(zip(list1, list2))
            f = lambda x:x.subs(subs_dict)
            ev = ev.applyfunc(f)
            REV = REV.applyfunc(f)
            LEV = LEV.applyfunc(f)
            ev_dict = {}
            REV_dict = {}
            LEV_dict = {}
            definitions ={}

            subs_list = [{Symbol('un_0'): 1, Symbol('un_1'): 0},
             {Symbol('un_0'): 0, Symbol('un_1'): 1}]
            equations = [Eq(a,b) for a,b in subs_dict.items()]
            eq_directions = {}
            for no,direction in enumerate(coordinates):
                g = lambda x:x.subs(subs_list[no]).simplify()
                definitions[direction] = ev.applyfunc(g).atoms(Symbol).union(REV.applyfunc(g).atoms(Symbol)).union(LEV.applyfunc(g).atoms(Symbol))
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)
                eq_directions[direction] = [e.subs(subs_list[no]) for e in equations]

        if ndim == 2:
            list1 = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            list2 = ['rho/(a*sqrt(2))', '1/(rho*a*sqrt(2))', 'k0*u0 + k1*u1', '(gama-1)*((u0**2 + u1**2)/2)', 'sqrt(k0**2 + k1**2)']
            list3 = ['k0_t','k1_t', 'theta_t', 'U']
            list4 = ['k0/k', 'k1/k', 'k0_t*u0 + k1_t*u1', 'k0*u0+k1*u1']

            list1 = [parse_expr(l,evaluate=False) for l in list1]
            list2 = [parse_expr(l, evaluate=False) for l in list2]
            list3 = [parse_expr(l, evaluate=False) for l in list3]
            list4 = [parse_expr(l, evaluate=False) for l in list4]
            subs_dict = dict(zip(list1, list2))
            subs_dict2 = dict(zip(list3, list4))

            ev = 'diag([U, U, U+a*k, U-a*k])'
            ev = parse_expr(ev, evaluate=False)
            REV = 'Matrix([[1,0,alpha,alpha], [u0,k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [u1,-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [phi_sq/(gama-1),rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + a*theta_t),alpha*((phi_sq+a**2)/(gama-1) - a*theta_t)]])'
            REV = parse_expr(REV, evaluate=False)
            LEV = 'Matrix([[1-phi_sq/a**2,(gama-1)*u0/a**2,(gama-1)*u1/a**2,-(gama-1)/a**2], [-(k1_t*u0-k0_t*u1)/rho,k1_t/rho,-k0_t/rho,0], [bta*(phi_sq-a*theta_t),bta*(k0_t*a-(gama-1)*u0),bta*(k1_t*a-(gama-1)*u1),bta*(gama-1)], [bta*(phi_sq+a*theta_t),-bta*(k0_t*a+(gama-1)*u0),-bta*(k1_t*a+(gama-1)*u1),bta*(gama-1)]])'
            LEV = parse_expr(LEV, evaluate=False)

            # Apply the sub
            f = lambda x:x.subs(subs_dict, evaluate=False)
            g = lambda x:x.subs(subs_dict2, evaluate=False)
            ev = ev.applyfunc(f).applyfunc(g).applyfunc(g)
            REV = REV.applyfunc(f).applyfunc(g).applyfunc(g)
            LEV = LEV.applyfunc(f).applyfunc(g).applyfunc(g)

            # Definitions required for each direction are
            definitions ={}

            # Dictionaries to store ev, REV and LEV for each direction
            ev_dict = {}
            REV_dict = {}
            LEV_dict = {}
            # CHANGE THIS to block dependent when block has the metric features
            metric_transformations = False
            if not metric_transformations:
                fact1 = 1
                fact2 = 1
            if metric_transformations:
                fact1 = ((met_symbols[0,0])**2 + (met_symbols[0,1])**2)**(Rational(1,2))
                fact2 = ((met_symbols[1,0])**2 + (met_symbols[1,1])**2)**(Rational(1,2))
            subs_list = [{Symbol('k0'): met_symbols[0,0], Symbol('k1'): met_symbols[0,1], Symbol('k'): fact1, Symbol('gama'): Rational(7,5)},
                        {Symbol('k0'): met_symbols[1,0], Symbol('k1'): met_symbols[1,1], Symbol('k'): fact2, Symbol('gama'): Rational(7,5)}]
            # equations = [Eq(a,b) for a,b in subs_dict.items()]
            eq_directions = {}
            for no,direction in enumerate(coordinates):
                g = lambda x:x.subs(subs_list[no], evaluate=False)
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)

        elif ndim == 3:
            list1 = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            list2 = ['rho/(a*sqrt(2))', '1/(rho*a*sqrt(2))', 'k0*u0 + k1*u1 + k2*u2', '(gama-1)*((u0**2 + u1**2 + u2**2)/2)', 'sqrt(k0**2 + k1**2 + k2**2)']
            list3 = ['k0_t','k1_t', 'k2_t', 'theta_t', 'U']
            list4 = ['k0/k', 'k1/k', 'k2/k', 'k0_t*u0 + k1_t*u1 + k2_t*u2', 'k0*u0+k1*u1+k2*u2']

            list1 = [parse_expr(l,evaluate=False) for l in list1]
            list2 = [parse_expr(l, evaluate=False) for l in list2]
            list3 = [parse_expr(l, evaluate=False) for l in list3]
            list4 = [parse_expr(l, evaluate=False) for l in list4]
            subs_dict = dict(zip(list1, list2))
            subs_dict2 = dict(zip(list3, list4))

            ev = 'diag([U, U, U, U+a*k, U-a*k])'
            ev = parse_expr(ev, evaluate=False)
            REV = 'Matrix([[k0_t,k1_t,k2_t,alpha,alpha], [k0_t*u0,k1_t*u0-k2_t*rho,k2_t*u0+k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [k0_t*u1+k2_t*rho,k1_t*u1,k2_t*u1-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [k0_t*u2-k1_t*rho,k1_t*u2+k0_t*rho,k2_t*u2,alpha*(u2+k2_t*a),alpha*(u2-k2_t*a)], [k0_t*phi_sq/(gama-1) + rho*(k2_t*u1-k1_t*u2),k1_t*phi_sq/(gama-1) + rho*(k0_t*u2-k2_t*u0),k2_t*phi_sq/(gama-1) + rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + theta_t*a),alpha*((phi_sq+a**2)/(gama-1) - theta_t*a)]])'
            REV = parse_expr(REV, evaluate=False)
            LEV = 'Matrix([[k0_t*(1-phi_sq/a**2)-(k2_t*u1-k1_t*u2)/rho,k0_t*(gama-1)*u0/a**2,k0_t*(gama-1)*u1/a**2 + k2_t/rho,k0_t*(gama-1)*u2/a**2 - k1_t/rho,-k0_t*(gama-1)/a**2], [k1_t*(1-phi_sq/a**2)-(k0_t*u2-k2_t*u0)/rho,k1_t*(gama-1)*u0/a**2 - k2_t/rho,k1_t*(gama-1)*u1/a**2,k1_t*(gama-1)*u2/a**2 + k0_t/rho,-k1_t*(gama-1)/a**2], [k2_t*(1-phi_sq/a**2) - (k1_t*u0-k0_t*u1)/rho,k2_t*(gama-1)*u0/a**2 + k1_t/rho,k2_t*(gama-1)*u1/a**2 - k0_t/rho,k2_t*(gama-1)*u2/a**2,-k2_t*(gama-1)/a**2], [bta*(phi_sq-theta_t*a),-bta*((gama-1)*u0-k0_t*a),-bta*((gama-1)*u1-k1_t*a),-bta*((gama-1)*u2 - k2_t*a),bta*(gama-1)], [bta*(phi_sq+theta_t*a),-bta*((gama-1)*u0+k0_t*a),-bta*((gama-1)*u1+k1_t*a),-bta*((gama-1)*u2+k2_t*a),bta*(gama-1)]])'
            LEV = parse_expr(LEV, evaluate=False)

            # Apply the sub
            f = lambda x:x.subs(subs_dict, evaluate=False)
            g = lambda x:x.subs(subs_dict2, evaluate=False)
            ev = ev.applyfunc(f).applyfunc(g).applyfunc(g)
            REV = REV.applyfunc(f).applyfunc(g).applyfunc(g)
            LEV = LEV.applyfunc(f).applyfunc(g).applyfunc(g)

            # Definitions required for each direction are
            definitions ={}

            # Dictionaries to store ev, REV and LEV for each direction
            ev_dict = {}
            REV_dict = {}
            LEV_dict = {}
            metric_transformations = False

            if not metric_transformations:
                fact1 = 1
                fact2 = 1
                fact3 = 1
            if metric_transformations:
                fact1= sqrt(met_symbols[0,0]**2 + met_symbols[0,1]**2 + met_symbols[0,2]**2)
                fact2= sqrt(met_symbols[1,0]**2 + met_symbols[1,1]**2 + met_symbols[1,2]**2)
                fact3= sqrt(met_symbols[2,0]**2 + met_symbols[2,1]**2 + met_symbols[2,2]**2)
            subs_list = [{Symbol('k0'): met_symbols[0,0], Symbol('k1'): met_symbols[0,1], Symbol('k2'): met_symbols[0,2], Symbol('k'): fact1, Symbol('gama'): Rational(7,5)},
                        {Symbol('k0'): met_symbols[1,0], Symbol('k1'): met_symbols[1,1], Symbol('k2'): met_symbols[1,2], Symbol('k'): fact2, Symbol('gama'): Rational(7,5)},
                        {Symbol('k0'): met_symbols[2,0], Symbol('k1'): met_symbols[2,1], Symbol('k2'): met_symbols[2,2], Symbol('k'): fact2, Symbol('gama'): Rational(7,5)}]
            # equations = [Eq(a,b) for a,b in subs_dict.items()]
            eq_directions = {}
            for no,direction in enumerate(coordinates):
                g = lambda x:x.subs(subs_list[no], evaluate=False)
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)


        return ev_dict, LEV_dict, REV_dict