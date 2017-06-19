#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy

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
from .opensbliobjects import DataSetBase, GridIndexedBase, ConstantObject

from sympy.core import Symbol,S, symbols
from .datatypes import *


class GridVariable(Symbol):

    """ Defines a new symbolic variable on the Grid.
        This should be used to create thread or process local variables in a kernel.
        dtype (datatype) and is_constant arguments are set to None and False initially.
        Attributes
        ==========
        ``dtype`` : data type of the variable, defaults to None
        ``is_constant`` : weather the variable is a constant or not.
        Notes
        ==========
        Grid variable cannot be a constant but used to fit in the abstraction.
        Later one can use this to define variables in Fortran subroutine for example
        For various dtypes used in OpenSBLI see ``datatypes.py``
        Parameter:
        ==========
        arg: string: variable: Name of the variable.
        :returns  variable defined on the grid.
        Examples:
        ==========
        >>> a = GridVariable("variable1")
        >>> srepr(variable1)
        GridVariable('variable1')
        >>> a.__dict__
        {dtype: None, is_constant: False}
        """

    def __new__(self, variable):
        self = Symbol.__new__(self, str(variable))
        self.is_constant = False
        self.dtype = None
        return self


class WorkDataSet():
    """ Base object for using work arrays contains differnt attributes to control the flow of work arrays
    in the opensbli framework. This should be used in conjunction with a Simulation block
    see ``SimulationBlock`` in ``block.py``. 
    Instantiated from Grid
    Some hard coded dependencies are work array name is always ``wk`` appended by the work index
    """
    def __init__(self):
        self.work_name = 'wk%d' # Work array name
        self.work_index = 0 # Index of the work array, this is update when ever a new work array is called
        self.stored_index = 0
        self.dtype = None # Place holder to save dtype not used currently
        return

    def work_array(self, name = None, location=None):
        """ Sets up a opensbli DataSet, indexed by the relative location.
        By default the range of the data set is that of the the block
        Parameters
        ==========
        :arg str name: The desired name of the work array (optional), defaults to ``wk``
        :arg location: List of relative indexing for the array (optional), defaults to grid intersection location
        :returns: Opensbli Dataset defined on the grid
        :rtype: opensbli.DataSet
        Notes
        ==========
        If multiple work arrays are required then the work_index should be incremented
        before calling this function. If required optionally the user can store the
        work index in their routine and reset the work index at the end of their routine
        if all the computations to be performed and the evaluated work arrays are not
        used any more. This helps in reducing the total memory foot print of the generated
        code.
        Example
        ==========
        >>> block = SimulationBlock(2, blocknumber=1)
        >>> b.work_array()
        wk0_B1[0,0]
        >>> b.increase_work_index
        >>> b.work_array()
        wk1_B1[0,0]
        >>> b.reset_work_index
        >>> b.work_array()
        wk0_B1[0,0]
        >>> b.increase_work_index
        >>> b.store_work_index
        >>> b.work_array()
        wk1_B1[0,0]
        >>> b.increase_work_index
        >>> b.work_array()
        wk2_B1[0,0]
        >>> b.increase_work_index
        >>> b.work_array()
        wk3_B1[0,0]
        >>> b.reset_work_to_stored
        >>> b.work_array()
        wk1_B1[0,0]
        """
        if not name:
            base = self.work_name%self.work_index
        else:
            base = name
        ret = DataSetBase(base)
        if location:
            if len(location) != self.ndim:
                raise ValueError("")
        else:
            # Default location to the current grid location
            location = [0 for i in range(self.ndim)]
        ret = ret[location]
        return ret
    @property
    def reset_work_index(self):
        """Resets the work index to zero. Used when we want to re-use work arrays
        in the discretisation
        """
        self.work_index = 0
        return
    @property
    def increase_work_index(self):
        """Increments the work array index by 1. This helps in setting up the next work array
        index
        """
        self.work_index = self.work_index +1
        return
    @property
    def store_work_index(self):
        """Stores the current work array index
        """
        self.stored_index = self.work_index
        return
    @property
    def reset_work_to_stored(self):
        """Resets the work array index to the index when the last ``store_work_index`` is called
        """
        self.work_index = self.stored_index
        return

class Grid(WorkDataSet):

    """ The numerical grid for a block and contains grid parameters"""

    def __init__(self):
        """ Initialise the gridobject and its parameters. This is instantiated from SimulationBlock
        see, ``block.py``
        """

        # Define the control parameters, currently not used
        self.define_control_parameters()

        # Instantiate WorkDataSet
        WorkDataSet.__init__(self)
        
        shape = symbols('block%dnp0:%d'%(self.blocknumber, self.ndim), integer=True)
        # Opensbli symbolic number of points, ConstantObjects
        self.shape = [ConstantObject("%s"%s, integer=True) for s in shape]
        # Sympy Idx type shape, lower and upper can be modified if required
        self.Idxed_shape = [Idx(Symbol('i%d'%dim, integer = True),(0, self.shape[dim])) for dim in range(self.ndim)]
        # For easier access ranges are created
        self.ranges = [[s.lower, s.upper] for s in self.Idxed_shape]
        # Grid spaciing in the number of dimensions, these are of type Constant Object
        self.deltas = [ConstantObject("Delta%dblock%d"%(dire,self.blocknumber)) for dire in range(self.ndim)]
        # Add the constants to ConstantsToDeclare
        from .kernel import ConstantsToDeclare as CTD
        for d in self.deltas:
            CTD.add_constant(d)
        for s in self.shape:
            CTD.add_constant(s, dtype = Int())
        # Name for the grid indices access (instead if i,j,k we use idx[0:ndim])
        g = GridIndexedBase('idx', self)
        self.grid_indexes = [g[i] for i in range(self.ndim)]
        return
    
    def define_control_parameters(self):
        """Not used, these should be used for further optimisations
        """
        # Parameters for optimizations / different schemes for convective and viscous fluxes
        self.store_derivatives = True # By default store derivatives is set to true
        self.derivatives_to_store = set()
        self.group_derivatives = False # This is another way of reducing data transfer.
        self.local_variables = True
        self.cns_type_code = False # This writes the code similat to traditional compressible NS
        self.sbli_rhs_discretisation = False # this is a pilot implementation of SBLI RHS
        return