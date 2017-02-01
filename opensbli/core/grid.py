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

from sympy.tensor import IndexedBase, Idx
from .opensbliobjects import DataSet

from sympy.core import Symbol,S, symbols

class DataType(object):
    pass

class Double(DataType):
    def __init__(self):
        return

class Float(DataType):
    def __init__(self):
        return

class UserDefined(DataType):
    """ User defined datatype this is either float or double depending on input"""
    def __init__(self):
        return

class GridVariable(Symbol):

    """ A symbolic variable defined on the Grid. """

    def __new__(self, variable):
        self = Symbol.__xnew__(self, variable)
        self.is_constant = False
        self.dtype = UserDefined()
        return self

#from .opensbliobjects import DataSer
class WorkDataSet():

    def __init__(self):
        self.work_name = 'wk%d'
        self.work_index = 0
        self.stored_index = 0
        self.dtype = UserDefined()
        return

    def work_array(self, name = None, location=None):
        """ Sets up a work array indexed by the Grid directions.
        No shape information will be provided, since the shape of the arrays might change based on the computations (including halos or excluding halos).

        :arg str name: The desired name of the work array.
        :returns: The Indexed object representing the work array defined on the Grid.
        :rtype: sympy.Indexed
        """
        if not name:
            base = self.work_name%self.work_index
        else:
            base = name
        # By default the range of the DataSet is the range of the grid
        if location:
            ret = DataSet(base, **{'location':location})
        else:
            ret = DataSet(base)
        ##print([[S.Zero, self.shape[i]] for i,val in enumerate(self.shape)])
        #ret.set_ranges([[S.Zero, self.shape[i]] for i,val in enumerate(self.shape)])
        return ret
    @property
    def reset_work_index(self):
        self.work_index = 0
        return
    @property
    def increase_work_index(self):
        self.work_index = self.work_index +1
        return
    @property
    def store_work_index(self):
        self.stored_index = self.work_index
        return
    @property
    def reset_work_to_stored(self):
        self.work_index = self.stored_index
        return
#Idx.lower.setter()

#from .bcs import BoundaryConditionTypes

class Grid(WorkDataSet):

    """ The numerical grid of solution points on which to discretise the equations. """

    def __init__(self):
        """ Initialise the grid of dimension ndim, and number of points nx0 x nx1 x nx2 (for the case of ndim = 3).

        :arg int ndim: The dimension of the grid.
        :arg dict grid_data: Optional user-defined grid parameters including the exact number of points and grid point spacing. If not provided, symbolic representations are used instead.
        :returns: None
        """

        # Define the control parameters
        self.define_control_parameters()

        # Instantiate WorkDataSet
        WorkDataSet.__init__(self)
        inds = symbols('shape_0:%d'%self.ndim, integer=True)
        self.shape = symbols('block%dnp_0:%d'%(self.blocknumber, self.ndim), integer=True)
        self.idx_shapes = [Idx(Symbol('i%d'%dim, integer = True),(0, self.shape[dim])) for dim in range(self.ndim)]
        self.ranges = [[s.lower, s.upper] for s in self.idx_shapes]
        #self.shape = [Idx(i, (upper[no])) for no,i in enumerate(inds)]
        #print(self.shape[0].args)
        #self.shape[0].args[1] = tuple([self.shape[0].lower, upper[1]])
        #print(self.shape[0].lower, self.shape[0].upper)
        return
    #def @upper.

    def define_control_parameters(self):
        # Parameters for optimizations / different schemes for convective and viscous fluxes
        self.store_derivatives = True # By default store derivatives is set to true
        self.derivatives_to_store = set()
        self.group_derivatives = False # This is another way of reducing data transfer.
        self.local_variables = True
        self.cns_type_code = False # This writes the code similat to traditional compressible NS
        self.sbli_rhs_discretisation = False # this is a pilot implementation of SBLI RHS
        return