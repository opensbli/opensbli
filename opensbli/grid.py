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
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>

from sympy import *

from .equations import EinsteinTerm


class Grid(object):

    """ The numerical grid of solution points on which to discretise the equations. """

    def __init__(self, ndim, grid_data=None):
        """ Initialise the grid of dimension ndim, and number of points nx0 x nx1 x nx2 (for the case of ndim = 3).

        :arg int ndim: The dimension of the grid.
        :arg dict grid_data: Optional user-defined grid parameters including the exact number of points and grid point spacing. If not provided, symbolic representations are used instead.
        :returns: None
        """

        # Number of grid points in each dimension.
        self.shape = tuple(symbols('nx0:%d' % ndim, integer=True))

        # Indices of the grid solution points in each dimension.
        self.indices = tuple(Symbol('i%d' % ind, integer=True) for ind, val in enumerate(self.shape))

        self.uniform = [True for ind, val in enumerate(self.shape)]

        # Define the grid point spacing term.
        di = EinsteinTerm('deltai_i')
        di.is_constant = True

        # Grid point spacing in each dimension.
        self.deltas = di.get_array(di.get_indexed(len(self.shape)))

        # Halo points. This will be populated when the stencil is created on the Grid.
        self.halos = []

        # FIXME: This works fine for now. But need a better idea.
        self.Idx = [Idx('idx[%d]' % ind) for ind, val in enumerate(self.shape)]

        # Use user-define grid data, if available and store in a dictionary to use it later.
        self.grid_data_dictionary = {}
        if grid_data:
            variables = [str(d) for d in self.deltas] + [str(s) for s in self.shape]
            values = [grid_data['delta'][i] for i in range(ndim)] + [grid_data['number_of_points'][i] for i in range(ndim)]
            self.grid_data_dictionary = dict(zip(variables, values))
            self.total_points = 1.0
            for r in range(ndim):
                self.total_points = self.total_points*grid_data['number_of_points'][r]

        # Coordinate indices in space (e.g. x0, x1, x2)
        coordinates = EinsteinTerm('x_i')
        coordinates.is_constant = True
        # This is essentially the indexed array of the "x_i" EinsteinTerm.
        self.coordinates = coordinates.get_array(coordinates.get_indexed(len(self.shape))).tolist()

        # Generate a mapping between the grid indices (e.g. i0, i1, i2) and the coordinate indices (e.g. x0, x1, x2).
        # This removes the dependancy of the coordinates in the spatial discretisation or Diagnostics or anywhere else.
        self.mapped_indices = dict(zip(tuple(self.coordinates)+self.indices, self.indices+self.indices))
        self.mapped_indices[EinsteinTerm('t')] = ''  # Also handle the time 't' as a special case.
        # NOTE: Indices that are mapped onto the grid should be populated implicitly.

        return

    def work_array(self, name):
        """ Sets up a work array indexed by the Grid directions.
        No shape information will be provided, since the shape of the arrays might change based on the computations (including halos or excluding halos).

        :arg str name: The desired name of the work array.
        :returns: The Indexed object representing the work array defined on the Grid.
        :rtype: sympy.Indexed
        """

        base = IndexedBase('%s' % name)
        base.is_grid = True
        base.is_constant = False
        return base[self.coordinates]

    def get_array_on_grid(self, array):
        """ Create a new IndexedBase object and set the attribute is_grid to True.
        Returns the new array with same name and same indices as the input array.
        The reason for the work-around is it is better to create an IndexedBase object and add attributes to it.

        :arg sympy.Indexed array: The array to be converted onto the Grid.
        :returns: The Indexed object representing the array defined on the Grid.
        :rtype: sympy.Indexed
        """
        base = array.base
        base.is_grid = True
        base.is_constant = False
        return base[array.indices]

    def indexed_by_grid(self, variable):
        """ Convert a variable/function or Indexed object to an Indexed object indexed by the Grid indices.

        :arg variable: The variable to convert to a Grid-based Indexed variable
        :arg grid: The numerical Grid of solution points.
        :returns: An Indexed variable, which is the same variable as the one provided, but is indexed by the Grid indices.
        :rtype: sympy.Indexed
        """

        if isinstance(variable, Indexed):
            base = IndexedBase('%s' % variable.base)
        elif isinstance(variable, Function):
            base = IndexedBase('%s' % variable.func)
        else:
            raise ValueError("Only functions or Indexed Objects are supported", variable)
        base.is_grid = True
        base.is_constant = False
        return base[self.indices]

    def grid_variable(self, name):
        """ Define a variable on the grid. This is not an Indexed variable but varies with the grid.
        Can be used as a local variable to be defined in a kernel.

        :arg str name: The name of the variable
        :returns: The grid variable.
        :rtype: opensbli.GridVariable
        """
        variable = GridVariable(str(name))
        return variable


class GridVariable(Symbol):

    """ A symbolic variable defined on the Grid. """

    def __new__(self, variable):
        self = Symbol.__xnew__(self, variable)
        return self
