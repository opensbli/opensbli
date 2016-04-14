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


class ReductionVariable(Symbol):

    def __new__(self, var):
        self = Symbol.__xnew__(self, var)
        return self


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
            indexes = [vin.indices for vin in allindexed if vin.base == v]
            if grid:
                v = self.set_grid_arrays(v, grid, indexes)
            self.inputs[v] = indexes
        for v in indexbase_outs:
            indexes = [vout.indices for vout in allindexed if vout.base == v]
            if grid:
                v = self.set_grid_arrays(v, grid, indexes)
            self.outputs[v] = indexes
        for v in indexbase_inouts:
            indexes = [vinout.indices for vinout in allindexed if vinout.base == v]
            if grid:
                v = self.set_grid_arrays(v, grid, indexes)
            self.inputoutput[v] = indexes

        idxs = flatten([list(e.rhs.atoms(Idx)) for e in self.equations])
        if idxs:
            self.has_Idx = True
        else:
            self.has_Idx = False

        self.reductions = flatten([list(e.rhs.atoms(ReductionVariable)) for e in self.equations])
        self.gridvariable = flatten([list(e.lhs.atoms(GridVariable)) for e in self.equations])
        if grid:
            self.constants = set(consts).difference(grid.mapped_indices.keys())
        else:
            self.constants = set(consts)

        return

    def set_grid_arrays(self, array, grid, indexes):
        """
        Sets the Indexed object attribute is_grid to True if all the indices of an indexed object
        are in mapped_indices dictionary of the grid

        """
        ets = [list(ind) for ind in indexes]
        ets = [list(et.atoms(Symbol)) for et in flatten(ets)]
        ets = (set(flatten(ets)))
        if all(index in grid.mapped_indices.keys() for index in ets):
            array.is_grid = True
        else:
            array.is_grid = False
        return array
