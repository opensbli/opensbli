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

from sympy import *
from .equations import EinsteinTerm


def decreasing_order(s1, s2):
    return cmp(len(s2.args), len(s1.args))


def increasing_order(s1, s2):
    return cmp(len(s1.args), len(s2.args))


def group_derivatives(derivatives):
    """ For each SymPy Derivative term, get the term that it is operating on and
    formula (function, Derivative(function)) pairs in a dictionary.

    :arg list derivatives: A list of all the Derivative terms to consider.
    :returns: A dictionary of the (function, Derivative(function)) pairs.
    :rtype: dict
    """

    derivative_dict = {}

    for derivative in derivatives:
        if derivative.args[0] in derivative_dict.keys():
            derivative_dict[derivative.args[0]] += [derivative]
        else:
            derivative_dict[derivative.args[0]] = [derivative]

    for key, value in derivative_dict.iteritems():
        if len(value) > 1:
            derivative_dict[key] = (sorted(value, cmp=increasing_order))

    return derivative_dict


def get_indexed_variables(equations):
    """ Return all the variables in the equations that are of type Indexed.

    :arg list equations: A list of SymPy equations to consider.
    :returns: A list of the variables in the equations that are Indexed, and also the count of all the specific terms in the equations (Indexed or not).
    :rtype: (list, int)
    """

    variables = []
    count = {}
    for eq in equations:
        pot = preorder_traversal(eq)
        for p in pot:
            if p in variables:
                pot.skip()
                count[p] = count[p]+1
                continue
            elif isinstance(p, Indexed):
                pot.skip()
                variables.append(p)
                count[p] = 1
            else:
                continue
    return variables, count


def substitute_work_arrays(ordered_evaluations, evaluations, equations):
    updated_equations = [eq for eq in equations]
    for equation_number, equation in enumerate(updated_equations):
        spatial_derivatives = [ev for ev in ordered_evaluations if isinstance(ev, Derivative)]
        formulas = [ev for ev in ordered_evaluations if isinstance(ev, Indexed)]
        spatial_derivatives = (sorted(spatial_derivatives, cmp=decreasing_order))
        for var in spatial_derivatives + formulas:
            new = evaluations[var].work
            updated_equations[equation_number] = updated_equations[equation_number].subs(var, new)
    return updated_equations


def update_work_arrays(ordered_evaluations, evaluations, work_array_name, work_array_index, grid):
    forms = [ev for ev in ordered_evaluations if isinstance(ev, Indexed)]
    for ev in forms:
        evaluations[ev].work = ev
    # Update the work arrays for the derivatives
    derivatives = [ev for ev in ordered_evaluations if isinstance(ev, Derivative)]
    for der in derivatives:
        wk = grid.work_array('%s%d' % (work_array_name, work_array_index))
        work_array_index += 1
        evaluations[der].work = wk
    return evaluations, work_array_index


def group_formulas(formulas, evals, known):
    """ This groups the formulas. """
    ranges = [evals[ev].evaluation_range for ev in formulas]
    subevals = flatten([evals[ev].subevals for ev in formulas])
    subeval_truth = [ev is None for ev in subevals]
    range_truth = [ranges[0][i] == val[i] for val in ranges for i in range(len(ranges[0]))]
    eqs = [Eq(ev, evals[ev].formula) for ev in formulas]
    range_dictionary = dict(zip(eqs, ranges))

    grouped_eq = []
    non_group = []
    for number, form in enumerate(formulas):
        if all(ev in known for ev in evals[form].requires):
            if all(range_truth) and all(subeval_truth):
                grouped_eq += [eqs[number]]
            else:
                non_group += [eqs[number]]
        else:
            non_group += [eqs[number]]
    return grouped_eq, non_group, range_dictionary


class SymbolicDerivative(object):

    """ The symbolic spatial derivatives of an arbitrary function 'F'
    on the numerical grid with the provided spatial scheme.

    For a wall boundary condition this will have a dependency on the grid range. """

    def __init__(self, spatial_scheme, grid):
        """ Initialise the spatial derivative, which gives the equations
        of spatial Derivatives for the various combinations of the spatial scheme and order of accuracy.

        :arg spatial_scheme: The spatial discretisation scheme to use.
        :arg grid: The numerical grid of solution points.
        :returns: None
        """
        self.derivative_direction = grid.indices
        self.index_mapping = grid.mapped_indices
        self.deltas = grid.deltas
        self.points = spatial_scheme.points
        return

    def get_derivative_formula(self, derivative):
        """ Return the formula for the derivative function. """
        order = len(derivative.args[1:])
        indices = list(derivative.args[1:])
        if order == 1 or len(set(indices)) == 1:
            wrt = indices[0]
            gridpoints = [wrt + i for i in self.points]
            formula = apply_finite_diff(order, gridpoints, [derivative.expr.subs({wrt: x}) for x in gridpoints], wrt)
            d1 = self.derivative_direction.index(self.index_mapping[derivative.args[1]])
            delta = self.deltas[d1]
            formula = formula*pow(delta, -order)
        elif order == 2:
            # Do the derivative of derivative
            raise NotImplementedError("Derivatives of order == 2 are not implemented.")
        else:
            raise NotImplementedError("Derivatives of order > 2 are not implemented.")
        return formula

    def get_derivative(self, derivative):
        """ Return a tuple to which the derivative formula exists in
        the already-evaluated derivatives.

        :arg derivative: The derivative you want to get the formula for.
        """
        order = len(derivative.args[1:])
        indices = []
        for arg in derivative.args[1:]:
            indices = indices + [self.derivative_direction.index(self.index_mapping[arg])]
        general_formula = []
        subevals = []
        requires = []
        if order == 1 or len(set(indices)) == 1:
            general_formula += [order, tuple(indices)]
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += list(derivative.args[0].atoms(Indexed))
        else:
            if len(derivative.args[0].atoms(Indexed)) > 1:
                subevals += [derivative.args[0]]
                requires += list(derivative.args[0].atoms(Indexed))
            else:
                subevals += [None]
                requires += [derivative.args[0]]
            general_formula += [order-1, tuple([indices[-1]])]
            requires += [Derivative(derivative.args[0], *derivative.args[1:-1])]

        return general_formula, subevals, requires


def get_used_formulas(formulas, equations):
    """ Return the formulas used in the equations. """
    variables, count = get_indexed_variables(equations)

    formulas = dict(zip([form.lhs for form in formulas], [form.rhs for form in formulas]))

    used_formulas = [Eq(var, formulas[var]) for var in variables if var in formulas.keys()]

    return used_formulas


def get_derivatives(equations):
    """ Return all the spatial Derivative terms in the equations.
    Any equations involving Derivative objects in terms of the time 't' are handled separately.

    :arg equations: A list of equations to search.
    :returns: All of the spatial Derivative objects and all of the temporal Derivative objects.
    """

    derivatives = []
    time_derivatives = []

    for eq in equations:
        pot = preorder_traversal(eq)

        for p in pot:
            if p in derivatives:
                pot.skip()
                continue
            elif isinstance(p, Derivative):
                if all(arg != EinsteinTerm('t') for arg in p.args):
                    pot.skip()
                    derivatives.append(p)
                else:
                    pot.skip()
                    time_derivatives.append(p)
            else:
                continue

    return derivatives, time_derivatives


def str_print(expr):
    val = str(expr)

    # Replacements
    indexed_objects = expr.atoms(Indexed)
    indexed_objects_replaced = [str(v.base) for v in indexed_objects]
    replacements = dict(zip([str(v) for v in indexed_objects], indexed_objects_replaced))
    replacements['Derivative'] = 'D'
    replacements[','] = ''
    for key, value in replacements.iteritems():
        val = val.replace(key, value)
    return val
