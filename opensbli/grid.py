from sympy import *

from .equations import EinsteinTerm

class Grid(object):
    
    """ The numerical grid of solution points on which to discretise the equations.
    If grid data dictionary is provided it will update with the exact values, else
    symbolic representations are generated
    """

    def __init__(self, ndim, grid_data=None):
        """ Initialise the grid of dimension ndim, and number of points nx0 x nx1 x nx2 (for the case of ndim = 3).
        
        :arg int ndim: The dimension of the grid.
        :arg dict grid_data: Optional user-defined grid parameters including the number of points and grid point spacing.
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
        
        # Use user-define grid data, if available.
        if grid_data:
            self.deltas = [grid_data['delta'][i] for i in range(ndim)]
            self.shape = tuple([grid_data['number_of_points'][i] for i in range(ndim)])
        
        # Halo points. This will be populated when the stencil is created on the Grid.
        self.halos = []
        
        # FIXME: This works fine for now. But need a better idea.
        self.Idx = [Idx('idx[%d]' % ind) for ind, val in enumerate(self.shape)]
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
