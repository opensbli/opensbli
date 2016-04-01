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
class ReductionVariable(Symbol):
    def __new__(self,var):
        self = Symbol.__xnew__(self, var)
        return self

class Kernel(object):

    """ A computational kernel which will be executed over all the grid points. """

    def __init__(self, equations, ranges, computation):
        """ Set up the kernel. This object will:

        1. Write the kernel's calling function based on the desired language.
            For this we require (for OPSC):
                a. Name of the kernel to call
                b. The block on which it should execute, and the dimensions of the block
                c. The ins, outs, and inouts variables, and their stencil of access and data type
                d. Indices of the array, if required

        2. Write the computational kernel, which requires writing kernel header and the computations
            For this we require (for OPSC):
                a. The kernel header with ins, outs, inouts and indices, specified along with the data type.
                b. To write the computations the Indexed objects are to be modified.
                   This modification requires OPS_ACC values which can also be populated and grid indices are replaced with 0's.

        All in all we require:
            1. Name
            2. Block (updated from the call to ops_write)
            3. ins, outs, inouts and so on
            4. Stencils of access
        """

        self.computation_type = computation
        self.ranges = ranges # Range of the indices of the points the kernel iterates over.
        self.name = None # None generates automatic kernel name
        if isinstance(equations, list):
            self.equations = equations
        else:
            self.equations = [equations]

        self.inputs = {}
        self.outputs = {}
        self.inputoutput = {}

        self.classify_grid_objects()

        return

    def classify_grid_objects(self):

        """ Classify the individual terms in the kernel's equation(s)
        as inputs, outputs, or inputoutputs (i.e. both an input and an output). """

        ins = []
        outs = []
        inouts = []
        consts = []

        for eq in self.equations:
            ins = ins + list(eq.rhs.atoms(Indexed))
            outs = outs + list(eq.lhs.atoms(Indexed))
            consts = consts + [et for et in list(eq.atoms(EinsteinTerm)) if et.is_constant]

        inouts = set(outs).intersection(set(ins))
        ins = set(ins).difference(inouts)
        outs = set(outs).difference(inouts)
        indexbase_ins = set([v.base for v in ins])
        indexbase_outs = set([v.base for v in outs])
        indexbase_inouts = set([v.base for v in inouts])
        for v in indexbase_ins:
            indexes = [vin.indices for vin in ins if vin.base==v]
            self.inputs[v] = indexes
        for v in indexbase_outs:
            indexes = [vout.indices for vout in outs if vout.base==v]
            self.outputs[v] = indexes
        for v in indexbase_inouts:
            indexes = [vinout.indices for vinout in inouts if vinout.base==v]
            self.inputoutput[v] = indexes

        idxs = flatten([list(e.rhs.atoms(Idx)) for e in self.equations])
        if idxs:
            self.has_Idx = True
        else:
            self.has_Idx = False

        self.reductions = flatten([list(e.rhs.atoms(ReductionVariable)) for e in self.equations])

        self.constants = set(consts)

        return
