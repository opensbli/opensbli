from sympy import *

from .equations import EinsteinTerm

class Grid(object):
    
    """ The numerical grid of solution points on which to discretise the equations.
    If grid data dictionary is provided it will update with the exact values, else
    symbolic representations are generated
    """

    def __init__(self, ndim, grid_data=None):
        self.shape = tuple(symbols('nx0:%d' % ndim, integer=True))
        # The indices of the grid solution points
        self.indices = [Symbol('i%d' % ind, integer = True) for ind, val in enumerate(self.shape)]
        self.uniform = [True for ind, val in enumerate(self.shape)]
        term = EinsteinTerm('deltai_i')
        term.is_constant = True
        self.deltas = term.get_array(term.get_indexed(len(self.shape)))
        if grid_data:
            self.deltas = [grid_data['delta'][i] for i in range(ndim)]
            self.shape = tuple([grid_data['number_of_points'][i] for i in range(ndim)])
        self.halos = []
        # FIXME: This works fine now. But need a better idea
        self.Idx = [Idx('idx[%d]' % ind) for ind, val in enumerate(self.shape)]
        return
        
    def work_array(self, name):
        """ No shape information will be provided; as the shape of the arrays might change based
        on the computations (including Halos or excluding halos) """
        
        base = IndexedBase('%s' % name)
        base.is_grid = True
        base.is_constant = False
        return base[self.indices]
