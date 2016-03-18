from sympy import *

from .equations import EinsteinTerm

class Grid(object):
    
    """ The numerical grid of solution points on which to discretise the equations. """

    def __init__(self, ndim):
        self.shape = tuple(symbols('nx0:%d' % ndim, integer=True))
        # The indices of the grid solution points
        self.indices = [Symbol('i%d' % ind, integer = True) for ind, val in enumerate(self.shape)]
        self.uniform = [True for ind, val in enumerate(self.shape)]
        term = EinsteinTerm('deltai_i')
        term.is_constant = True
        self.deltas = term.get_array(term.get_indexed(len(self.shape)))
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
