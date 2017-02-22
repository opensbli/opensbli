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

# Import local utility functions
from sympy.tensor.array import MutableDenseNDimArray, tensorcontraction
from opensbli.core.opensbliobjects import *
from sympy import *
from opensbli.core.opensbliobjects import ConstantObject, EinsteinTerm
from opensbli.core.opensblifunctions import CentralDerivative
from opensbli.core.opensbliequations import OpenSBLIExpression
from sympy.parsing.sympy_parser import parse_expr
from opensbli.core.opensblifunctions import *

class EulerEquations(object):
    def __init__(self, ndim, weno=True, **kwargs):
        """ """
        self.weno = weno
        self.constants = []
        self.general_equations()
        self.formulas = []
        self.velocity_components()
        self.pressure_equation()
        self.speed_of_sound()
        self.ndim = ndim
        coordinate_symbol = "x"

        self.time_derivative_symbol = CoordinateObject('t')
        cart = CoordinateObject('%s_i'%(coordinate_symbol))
        self.space_derivative_symbols = [cart.apply_index(cart.indices[0], dim) for dim in range(self.ndim)]
        # self.eq_to_vector_form()
        # self.formula_to_dict()
        self.REV = {}
        self.LEV = {}
        self.EV= {}
        # self.conser_to_primitive()
        return

    def get_time_derivative(self, eqns):
        """
        Get the time derivatives to add to the vector notation dictionary.
        """
        time_deriv = []
        for deriv in eqns.atoms(TemporalDerivative):
            time_deriv.append(deriv)
        return time_deriv

    def group_by_direction(self, eqn):
        """
        Group the derivatives by direction, one equation at a time given to this function.
        """
        all_WDS = []
        all_WDS += eqn.atoms(WenoDerivative)
        all_WDS = list(set(all_WDS))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped

    def general_equations(self):
        '''General equations for the compressible Navier Stokes are written in here, depending on the
        equations required these will be modified at the initialisation of the class
        '''
        # Shock capturing enabled
        if self.weno == True:
            scheme = "**{\'scheme\':\'Weno\'}"
        else:
            scheme = "**{\'scheme\':\'Central\'}"
        mass = "Eq(Der(rho,t), - Conservative(rhou_j,x_j, %s))" % scheme
        momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j+KD(_i,_j)*p,x_j, %s))" % scheme
        energy = "Eq(Der(rhoE,t), - Conservative((rhoE+p)*u_j,x_j, %s) )" % scheme
        self.equations = [mass, momentum, energy]
        return

    def velocity_components(self):
        velocity = "Eq(u_i, rhou_i/rho)"
        self.formulas = self.formulas + [velocity]
        return

    def pressure_equation(self, conservative= True):
        # if not conservative:
        #     pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*rho*(u_j*u_j)))"
        # pressure = "Eq(p, (gama-1)*(rhoE - (1/2)*(rhou_j*rhou_j/rho)))"
        #WARNING: Added pressure with KD instead of two dummy as before, is this formula correct?
        pressure = "Eq(p, (gama-1)*(rhoE - rho*(1/2)*(KD(_i,_j)*u_i*u_j)))"
        self.constants += ['gama']
        self.formulas = self.formulas + [pressure]
        return

    def speed_of_sound(self,conservative= True):
        # asq = "Eq(asq, gama*(gama-1)*(rhoE/rho - (1/2)*(rhou_j*rhou_j/(rho*rho))))"
        asq = "Eq(asq, gama*(gama-1)*(rhoE/rho - (1/2)*(KD(_i,_j)*u_i*u_j)))"
        a = "Eq(a, sqrt(gama*(gama-1)*(rhoE/rho - (1/2)*(KD(_i,_j)*u_i*u_j))))"
        self.formulas = self.formulas + [asq]
        self.formulas = self.formulas + [a]
        return

    def eq_to_vector_form(self, eqns):
        # Flatten the system of equations
        flat_eqns = flatten(eqns)
        # Generate dictionary to group the terms in vector form for each coordinate
        self.vector_notation = {key:[] for key in self.space_derivative_symbols}
        t = self.time_derivative_symbol
        self.vector_notation[t] = []
        # Loop over the equations in the system
        for eq in flat_eqns:
            # Group terms in this equation based on the direction the differentiation is applied
            grouped_dictionary = self.group_by_direction(eq)
            # Loop over the directions
            for key, derivatives in grouped_dictionary.iteritems():
                coordinate = self.space_derivative_symbols[key]
                # Check if the coordinate exists in the vector notation dictionary
                if coordinate in self.vector_notation.keys():
                    # Add the function being differentiated to its position in the vector notation
                    self.vector_notation[coordinate] = self.vector_notation[coordinate] + [deriv.args[0] for deriv in derivatives]
            # Find time derivatives and add them into the vector expression
            grouped_dictionary[t] = self.get_time_derivative(eq)
            self.vector_notation[t] = self.vector_notation[t] + [deriv.args[0] for deriv in grouped_dictionary[t]]
        # Convert to matrices
        for key, eq in self.vector_notation.iteritems():
            self.vector_notation[key] = Matrix(eq)
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
        local_dict = { 'Symbol': EinsteinTerm, 'gama': ConstantObject('gama')}

        if ndim == 1:
            matrix_symbols = ['H']
            matrix_formulae = ['a**2/(gama-1) + u0**2/2']
            matrix_symbols = [parse_expr(l, local_dict=local_dict) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict) for l in matrix_formulae]


            ev = 'diag([u0-a, u0, u0+a])'
            REV = 'Matrix([[1,1,1], [u0-a,u0,u0+a], [H-u0*a,u0**2 /2,H+u0*a]])'
            LEV = 'Matrix([[ u0*(2*H + a*u0 - u0**2)/(2*a*(2*H - u0**2)), (-H - a*u0 + u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)],[2*(H - u0**2)/(2*H - u0**2),2*u0/(2*H - u0**2), 2/(-2*H + u0**2)],[u0*(-2*H + a*u0 + u0**2)/(2*a*(2*H - u0**2)),  (H - a*u0 - u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)]])'
            ev = parse_expr(ev, local_dict=local_dict)
            REV = parse_expr(REV, local_dict=local_dict)
            LEV = parse_expr(LEV, local_dict=local_dict)


            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            f = lambda x:x.subs(subs_dict)
            ev = ev.applyfunc(f)
            REV = REV.applyfunc(f)
            LEV = LEV.applyfunc(f)
            ev_dict = {}
            REV_dict = {}
            LEV_dict = {}
            definitions ={}

            subs_list = [{EinsteinTerm('un_0'): 1, EinsteinTerm('un_1'): 0},
             {EinsteinTerm('un_0'): 0, EinsteinTerm('un_1'): 1}]
            equations = [Eq(a,b) for a,b in subs_dict.items()]
            eq_directions = {}
            for no,coordinate in enumerate(coordinates):
                direction =coordinate.direction
                g = lambda x:x.subs(subs_list[no]).simplify()
                definitions[direction] = ev.applyfunc(g).atoms(EinsteinTerm).union(REV.applyfunc(g).atoms(EinsteinTerm)).union(LEV.applyfunc(g).atoms(EinsteinTerm))
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)
                eq_directions[direction] = [e.subs(subs_list[no]) for e in equations]

        if ndim == 2:
            matrix_symbols = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            matrix_formulae = ['rho/(a*sqrt(2))', '1/(rho*a*sqrt(2))', 'k0*u0 + k1*u1', '(gama-1)*((u0**2 + u1**2)/2)', 'sqrt(k0**2 + k1**2)']
            matrix_symbols_2 = ['k0_t','k1_t', 'theta_t', 'U']
            matrix_formulae_2 = ['k0/k', 'k1/k', 'k0_t*u0 + k1_t*u1', 'k0*u0+k1*u1']

            matrix_symbols = [parse_expr(l, local_dict=local_dict) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict) for l in matrix_formulae]
            matrix_symbols_2 = [parse_expr(l, local_dict=local_dict) for l in matrix_symbols_2]
            matrix_formulae_2 = [parse_expr(l, local_dict=local_dict) for l in matrix_formulae_2]
            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            subs_dict2 = dict(zip(matrix_symbols_2, matrix_formulae_2))

            ev = 'diag([U, U, U+a*k, U-a*k])'
            REV = 'Matrix([[1,0,alpha,alpha], [u0,k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [u1,-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [phi_sq/(gama-1),rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + a*theta_t),alpha*((phi_sq+a**2)/(gama-1) - a*theta_t)]])'
            LEV = 'Matrix([[1-phi_sq/a**2,(gama-1)*u0/a**2,(gama-1)*u1/a**2,-(gama-1)/a**2], [-(k1_t*u0-k0_t*u1)/rho,k1_t/rho,-k0_t/rho,0], [bta*(phi_sq-a*theta_t),bta*(k0_t*a-(gama-1)*u0),bta*(k1_t*a-(gama-1)*u1),bta*(gama-1)], [bta*(phi_sq+a*theta_t),-bta*(k0_t*a+(gama-1)*u0),-bta*(k1_t*a+(gama-1)*u1),bta*(gama-1)]])'
            ev = parse_expr(ev, local_dict=local_dict)
            REV = parse_expr(REV, local_dict=local_dict)
            LEV = parse_expr(LEV, local_dict=local_dict)

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
            subs_list = [{EinsteinTerm('k0'): met_symbols[0,0], EinsteinTerm('k1'): met_symbols[0,1], EinsteinTerm('k'): fact1},
                        {EinsteinTerm('k0'): met_symbols[1,0], EinsteinTerm('k1'): met_symbols[1,1], EinsteinTerm('k'): fact2}]
            # equations = [Eq(a,b) for a,b in subs_dict.items()]
            eq_directions = {}
            for no,coordinate in enumerate(coordinates):
                direction =coordinate.direction
                g = lambda x:x.subs(subs_list[no], evaluate=False)
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)
            
        elif ndim == 3:
            matrix_symbols = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            matrix_formulae = ['rho/(a*sqrt(2))', '1/(rho*a*sqrt(2))', 'k0*u0 + k1*u1 + k2*u2', '(gama-1)*((u0**2 + u1**2 + u2**2)/2)', 'sqrt(k0**2 + k1**2 + k2**2)']
            matrix_symbols_2 = ['k0_t','k1_t', 'k2_t', 'theta_t', 'U']
            matrix_formulae_2 = ['k0/k', 'k1/k', 'k2/k', 'k0_t*u0 + k1_t*u1 + k2_t*u2', 'k0*u0+k1*u1+k2*u2']

            matrix_symbols = [parse_expr(l, local_dict=local_dict) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict) for l in matrix_formulae]
            matrix_symbols_2 = [parse_expr(l, local_dict=local_dict) for l in matrix_symbols_2]
            matrix_formulae_2 = [parse_expr(l, local_dict=local_dict) for l in matrix_formulae_2]
            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            subs_dict2 = dict(zip(matrix_symbols_2, matrix_formulae_2))

            ev = 'diag([U, U, U, U+a*k, U-a*k])'
            REV = 'Matrix([[k0_t,k1_t,k2_t,alpha,alpha], [k0_t*u0,k1_t*u0-k2_t*rho,k2_t*u0+k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [k0_t*u1+k2_t*rho,k1_t*u1,k2_t*u1-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [k0_t*u2-k1_t*rho,k1_t*u2+k0_t*rho,k2_t*u2,alpha*(u2+k2_t*a),alpha*(u2-k2_t*a)], [k0_t*phi_sq/(gama-1) + rho*(k2_t*u1-k1_t*u2),k1_t*phi_sq/(gama-1) + rho*(k0_t*u2-k2_t*u0),k2_t*phi_sq/(gama-1) + rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + theta_t*a),alpha*((phi_sq+a**2)/(gama-1) - theta_t*a)]])'
            LEV = 'Matrix([[k0_t*(1-phi_sq/a**2)-(k2_t*u1-k1_t*u2)/rho,k0_t*(gama-1)*u0/a**2,k0_t*(gama-1)*u1/a**2 + k2_t/rho,k0_t*(gama-1)*u2/a**2 - k1_t/rho,-k0_t*(gama-1)/a**2], [k1_t*(1-phi_sq/a**2)-(k0_t*u2-k2_t*u0)/rho,k1_t*(gama-1)*u0/a**2 - k2_t/rho,k1_t*(gama-1)*u1/a**2,k1_t*(gama-1)*u2/a**2 + k0_t/rho,-k1_t*(gama-1)/a**2], [k2_t*(1-phi_sq/a**2) - (k1_t*u0-k0_t*u1)/rho,k2_t*(gama-1)*u0/a**2 + k1_t/rho,k2_t*(gama-1)*u1/a**2 - k0_t/rho,k2_t*(gama-1)*u2/a**2,-k2_t*(gama-1)/a**2], [bta*(phi_sq-theta_t*a),-bta*((gama-1)*u0-k0_t*a),-bta*((gama-1)*u1-k1_t*a),-bta*((gama-1)*u2 - k2_t*a),bta*(gama-1)], [bta*(phi_sq+theta_t*a),-bta*((gama-1)*u0+k0_t*a),-bta*((gama-1)*u1+k1_t*a),-bta*((gama-1)*u2+k2_t*a),bta*(gama-1)]])'
            ev = parse_expr(ev, local_dict=local_dict)
            REV = parse_expr(REV, local_dict=local_dict)
            LEV = parse_expr(LEV, local_dict=local_dict)

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
            for no,coordinate in enumerate(coordinates):
                direction = coordinate.direction
                g = lambda x:x.subs(subs_list[no], evaluate=False)
                ev_dict[direction] = diag(*list(ev.applyfunc(g)))
                REV_dict[direction] = REV.applyfunc(g)
                LEV_dict[direction] = LEV.applyfunc(g)
        return ev_dict, LEV_dict, REV_dict