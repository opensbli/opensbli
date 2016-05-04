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
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>

from sympy import *

from .equations import EinsteinTerm
from .grid import GridVariable
from .utils import *


class Kernel(object):

    """ A computational kernel which will be executed over all the grid points. """

    def __init__(self, equations, ranges, computation, grid=None):
        """ Set up the computational kernel"""

        self.computation_type = computation
        self.ranges = ranges  # Range of the indices of the points the kernel iterates over.
        self.name = None  # None generates automatic kernel name
        if isinstance(equations, list):
            self.equations = equations
        else:
            self.equations = [equations]

        self.inputs = {}
        self.outputs = {}
        self.inputoutput = {}

        self.classify_grid_objects(grid)

        return

    def classify_grid_objects(self, grid):
        """ Classify the individual terms in the kernel's equation(s)
        as inputs, outputs, or inputoutputs (i.e. both an input and an output). """

        ins = []
        outs = []
        consts = []
        allindexed = []

        for eq in self.equations:
            ins = ins + list(eq.rhs.atoms(IndexedBase))
            outs = outs + list(eq.lhs.atoms(IndexedBase))
            allindexed = allindexed + list(eq.atoms(Indexed))
            consts = consts + [et for et in list(eq.atoms(EinsteinTerm)) if et.is_constant]

        indexbase_inouts = set(outs).intersection(set(ins))
        indexbase_ins = set(ins).difference(indexbase_inouts)
        indexbase_outs = set(outs).difference(indexbase_inouts)

        for v in indexbase_ins:
            indexes = list(set([vin.indices for vin in allindexed if vin.base == v]))
            if grid:
                v = self.set_grid_arrays(v, grid, indexes)
            self.inputs[v] = indexes
        for v in indexbase_outs:
            indexes = list(set([vout.indices for vout in allindexed if vout.base == v]))
            if grid:
                v = self.set_grid_arrays(v, grid, indexes)
            self.outputs[v] = indexes
        for v in indexbase_inouts:
            indexes = list(set([vinout.indices for vinout in allindexed if vinout.base == v]))
            if grid:
                v = self.set_grid_arrays(v, grid, indexes)
            self.inputoutput[v] = indexes

        idxs = flatten([list(e.rhs.atoms(Idx)) for e in self.equations])
        if idxs:
            self.has_Idx = True
        else:
            self.has_Idx = False

        from .diagnostics import ReductionVariable
        self.reductions = flatten([list(e.rhs.atoms(ReductionVariable)) for e in self.equations])
        self.gridvariable = flatten([list(e.lhs.atoms(GridVariable)) for e in self.equations])
        if grid:
            self.constants = set(consts).difference(grid.mapped_indices.keys())
        else:
            self.constants = set(consts)

        return

    def set_grid_arrays(self, array, grid, indexes):
        """ Sets the Indexed object attribute is_grid to True if all the indices of an Indexed object
        are in the 'mapped_indices' dictionary of the Grid. """
        ets = [list(ind) for ind in indexes]
        ets = [list(et.atoms(Symbol)) for et in flatten(ets)]
        ets = (set(flatten(ets)))
        if all(index in grid.mapped_indices.keys() for index in ets):
            array.is_grid = True
        else:
            array.is_grid = False
        return array


def create_derivative_kernels(derivatives, evals, spatial_derivative, work_array_name, work_array_index, grid):
    computations = []
    ranges = [evals[ev].evaluation_range for ev in derivatives]
    subevals = [evals[ev].subevals for ev in derivatives]
    require = [evals[ev].requires for ev in derivatives]
    for number, derivative in enumerate(derivatives):
        if not any(isinstance(req, Derivative) for req in require[number]):
            if all(subev is None for subev in subevals[number]):
                rhs = spatial_derivative.get_derivative_formula(derivative)
                eq = Eq(evals[derivative].work, rhs)
                name = str_print(derivative)
                computations.append(Kernel(eq, ranges[number], name, grid))
            else:
                # Store into temporary array the sub evaluation
                eqs = []
                temp_work_array_index = work_array_index
                for subev in subevals[number]:
                    wk = grid.work_array('%s%d' % (work_array_name, temp_work_array_index))
                    temp_work_array_index += 1
                    for req in require[number]:
                        local_range = evals[req].evaluation_range
                        subev = subev.subs(req, evals[req].work)
                    eqs.append(Eq(wk, subev))
                name = str_print(subev)
                computations.append(Kernel(eqs, local_range, name, grid))
                for eq in eqs:
                    new_derivative = derivative.subs(eq.rhs, eq.lhs)
                rhs = spatial_derivative.get_derivative_formula(new_derivative)
                eq = Eq(evals[derivative].work, rhs)
                name = str_print(derivative)
                computations.append(Kernel(eq, ranges[number], name, grid))
        else:
            new_derivative = derivative
            if all(subev is None for subev in subevals[number]):
                for req in require[number]:
                    new_derivative = new_derivative.subs(req, evals[req].work)
            else:
                raise NotImplementedError("Sub-evaluations in a mixed derivative", grid)
            rhs = spatial_derivative.get_derivative_formula(new_derivative)
            eq = Eq(evals[derivative].work, rhs)
            name = str_print(derivative)
            computations.append(Kernel(eq, ranges[number], name, grid))
    return computations


def create_formula_kernels(ordered_evaluations, evaluations, known, grid):
    computation_kernels = []
    forms = [ev for ev in ordered_evaluations if isinstance(ev, Indexed) and ev not in known]
    grouped, non_group, range_dictionary = group_formulas(forms, evaluations, known)
    if grouped:
        computation_kernels += [Kernel(grouped, range_dictionary[grouped[0]], "Grouped Formula Evaluation", grid)]
    for eq in non_group:
        computation_kernels += [Kernel(eq, range_dictionary[eq], "Non-Grouped Formula Evaluation", grid)]
    return computation_kernels
