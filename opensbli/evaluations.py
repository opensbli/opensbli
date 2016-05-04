#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

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
from .utils import *


class Evaluations(object):

    """ The evaluation of a LHS and RHS, containing information about what the LHS and RHS requires, whether there are any
    subevaluations that need doing (e.g. for second derivatives, e.g. d(du/dx)dy, du/dx should be evaluated first),
    and what the work arrays are. """

    def __init__(self, lhs, rhs, requires, subevals=None, wk=None):
        """ Set up the evaluation information. """

        if isinstance(lhs, Derivative):
            self.is_derivative = True
            self.is_formula = False

            if subevals:
                self.subevals = subevals
            else:
                self.subevals = [None]
            if wk:
                self.work = wk
            else:
                self.work = None

            self.formula = rhs
            self.requires = requires
            self.evaluation_range = []
        else:
            self.is_formula = True
            self.is_derivative = False

            if subevals:
                self.subevals = subevals
            else:
                self.subevals = [None]
            if wk:
                self.work = wk
            else:
                self.work = None

            self.formula = rhs
            self.requires = requires
            self.evaluation_range = []

        return


def create_derivative_evaluations(spatial_derivatives, evals, symbolic_derivative):
    """
    Derivative computations are evaluated seperately as they sometimes require evaluation of
    temporary work arrays
    """
    for out in spatial_derivatives:
        general_formula, subevals, requires = symbolic_derivative.get_derivative(out)
        evaluated = Evaluations(out, general_formula, requires, subevals, None)
        evals[out] = evaluated

    return evals


def create_formula_evaluations(all_formulas, evals):
    """
    Creates computation class for all the equations provided.
    This is used for creating any computation (Diagnostics, Spatial, Temporal, Shock-capturing)
    """
    for formula in all_formulas:
        evaluated = Evaluations(formula.lhs, formula.rhs, list(formula.rhs.atoms(Indexed)), None, None)
        evals[formula.lhs] = evaluated
    return evals


def sort_evaluations(order, evaluations, typef):
    """ Sort the evaluations based on the requirements of each term. For example, if we have
    the primitive variables p, u0, u1, and T, then the pressure p may depend on the velocity u0 and u1, and T may depend on p,
    so we need this be evaluate in the following order: u0, u1, p, T.


    :arg list order: The list of terms to sort.
    :arg evaluations: The evaluation information, containing dependency information.
    :arg typef: The type of term to sort.
    :returns: A list of ordered terms.
    :rtype: list
    """
    key_list = [key for key in evaluations.keys() if isinstance(key, typef) and key not in order]
    requires_list = [evaluations[key].requires for key in key_list]
    zipped = zip(key_list, requires_list)
    # Breaks after 1000 iterations
    iter_count = 0
    while key_list:
        iter_count = iter_count+1
        order += [x for (x, y) in zipped if all(req in order for req in y)]
        key_list = [key for key in evaluations.keys() if isinstance(key, typef) and key not in order]
        requires_list = [evaluations[key].requires for key in key_list]
        zipped = zip(key_list, requires_list)
        if iter_count > 1000:
            pprint([req for req in requires_list[0]])
            pprint(order)
            pprint([req in order for req in requires_list[0]])
            raise ValueError("Exiting sort evaluations ")

    return order


def set_range_of_evaluations(order_of_evaluations, evaluations, grid):
    """ Set the evaluation ranges of each Evaluation object based on the shape of the grid of points.
    First the ranges of derivatives are updated, then other ranges are updated. """

    derivatives = []
    for ev in order_of_evaluations:
        if isinstance(ev, Derivative):
            derivatives.append(ev)
        evaluations[ev].evaluation_range = [tuple([0, s]) for s in grid.shape]

    # Update the range for the derivatives
    grouped_derivatives = group_derivatives(derivatives)
    for key, value in grouped_derivatives.iteritems():
        for val in value:
            require = evaluations[val].requires
            formula = evaluations[val].formula
            direction = formula[1][0]
            halos = grid.halos[direction]
            for req in require:
                erange = list(evaluations[req].evaluation_range[direction])
                if erange[0] == 0 and erange[1] == grid.shape[direction]:
                    erange[0] += halos[0]
                    erange[1] += halos[1]
                evaluations[req].evaluation_range[direction] = tuple(erange)

    # Update the range for the formulas this might require some modifications
    for ev in order_of_evaluations:
        if isinstance(ev, Indexed):
            require = evaluations[ev].requires
            if require:
                for req in require:
                    evaluations[req].evaluation_range = evaluations[ev].evaluation_range

    return
