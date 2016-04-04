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
        self.indices = tuple(Symbol('i%d' % ind, integer = True) for ind, val in enumerate(self.shape))
        
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
            self.grid_data_dictionary  = dict(zip([str(d) for d in self.deltas], \
                [grid_data['delta'][i] for i in range(ndim)])) 
            self.shape = tuple([grid_data['number_of_points'][i] for i in range(ndim)])
        # Also require a mapping between grid indices and the coordinate directions
        # This is like the indexed array of x_i EinsteinTerm, this removes the dependancy of coordinate 
        # in the spatial descritisation or Diagnostics or any where else later this should be input to this
        
        di = EinsteinTerm('x_i')
        di.is_constant = True
        
        # Grid point spacing in each dimension.
        self.mapedindices = di.get_array(di.get_indexed(len(self.shape))).tolist()
        self.mapped_indices = dict(zip(self.mapedindices, self.indices))
        self.mapped_indices[EinsteinTerm('t')] = ''
        # indices that are mapped on to the grid should be populated Implicitly 
        return
        
    def work_array(self, name):
        """ Sets up a work array indexed by the Grid.
        No shape information will be provided, since the shape of the arrays might change based on the computations (including halos or excluding halos).
        
        :arg str name: The desired name of the work array.
        :returns: The Indexed object representing the work array defined on the Grid.
        :rtype: sympy.Indexed
        """
        
        base = IndexedBase('%s' % name)
        base.is_grid = True
        base.is_constant = False
        return base[self.indices]
