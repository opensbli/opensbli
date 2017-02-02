#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy and others

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

# @author: New structure implemented by Satya P Jammy (October, 2016)

from .grid import *
from sympy.matrices import *
from .bcs import BoundaryConditionTypes
from .opensbliequations import SimulationEquations

class KernelCounter():
    # Counter for the kernels
    def __init__(self):
        self.kernel_counter = 0

    def reset_kernel_counter(self):
        self.kernel_counter = 0
        return
    @property
    def increase_kernel_counter(self):
        self.kernel_counter = self.kernel_counter +1
        return

    def store_kernel_counter(self):
        self.stored_counter = self.kernel_counter
        return

    def reset_kernel_to_stored(self):
        self.kernel_counter = self.stored_counter
        return

class RationalCounter():
    # Counter for the kernels
    def __init__(self):
        self.name = 'rc%d' 
        self.rational_counter = 0

    @property
    def increase_rational_counter(self):
        self.rational_counter = self.rational_counter +1
        return
    @property
    def get_next_rational_constant(self):
        name = self.name % self.rational_counter
        self.increase_rational_counter
        ret = ConstantObject(name)
        return ret
class SimulationBlock(Grid, KernelCounter, BoundaryConditionTypes): # BoundaryConditionTypes add this later
    def __init__(self, ndim, block_number = None):
        if block_number:
            self.blocknumber = block_number
        else:
            self.blocknumber = 0
        self.ndim = ndim
        KernelCounter.__init__(self)
        Grid.__init__(self)
        RationalCounter.__init__(self)
        self.boundary_halos = [[set(), set()] for d in range(self.ndim)]
        self.block_datasets = {}
        self.constants = set()
        self.Rational_constants = {}
        return

    def set_block_number(self, number):
        self.blocknumber = number
        return

    def set_block_boundaries(self, bclist):
        self.set_boundary_types(bclist)

    def set_block_boundary_halos(self, direction, side, types):
        self.boundary_halos[direction][side].add(types)
        return

    def discretise(self):
        """
        In this the discretisation of the schemes in the list of equations is applied
        :arg list_of_equations: a list of the type of equations (simulation equations, Constituent relations,
        Metric equations, diagnostic equations etc)
        : ar
        """
        # perform the spatial discretisation of the equations using schemes
        for eq in self.list_of_equation_classes:
            eq.spatial_discretisation(self.discretisation_schemes, self)
            eq.apply_boundary_conditions(self)
        # Get the classes for the constituent relations
        crs = self.get_constituent_equation_class
        for clas in self.list_of_equation_classes:
            if clas not in crs:
                print clas, "Exitting"
                #exit()
        # perform the temporal discretisation of the equations for all equation classes
        # Later move TD to equations.td
        temporal = self.get_temporal_schemes
        for t in temporal:
            for eq in self.list_of_equation_classes:
                self.discretisation_schemes[t.name].discretise(eq, self)
        return

    def apply_boundary_conditions(self, arrays):
        kernels = []
        for no,b in enumerate(self.boundary_types):
            kernels += [self.apply_bc_direction(no, 0, arrays)]
            kernels += [self.apply_bc_direction(no, 1, arrays)]
        return kernels

    def apply_bc_direction(self, direction, side, arrays):
        kernel = self.boundary_types[direction][side].apply(arrays, direction, side, self)
        return kernel

    def set_equations(self, list_of_equations):
        self.list_of_equation_classes = list_of_equations
        return

    def set_discretisation_schemes(self, schemes):
        self.discretisation_schemes = schemes
        return

    @property
    def get_constituent_equation_class(self):
        from .opensbliequations import ConstituentRelations as CR
        CR_classes = []
        for sc in self.list_of_equation_classes:
            if isinstance(sc, CR):
                CR_classes += [sc]
        return CR_classes

    @property
    def get_temporal_schemes(self):
        temporal = []
        for sc in  self.discretisation_schemes:
            if self.discretisation_schemes[sc].schemetype == "Temporal":
                temporal += [self.discretisation_schemes[sc]]
        return temporal
    @property
    def collect_all_spatial_kernels(self):
        all_kernels = []
        for scheme in self.get_temporal_schemes:
            for key, value in scheme.solution.iteritems(): # These are equation classes
                if key.order >=0 and key.order <100: #Checks if the equation classes are part of the time loop
                    all_kernels += key.all_spatial_kernels
                else:
                    print 'NOPE' # Just checking
        return all_kernels

    def grid_generation(self):

        return

    def initial_conditions(self):

        return

    def io(self):
        return

    def pre_process_eq(self, eq_class):
        """
        These are type non Simulation equations
        """
        return

    def post_process_eq(self, eq_class_list):

        return
