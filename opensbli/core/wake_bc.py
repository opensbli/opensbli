from sympy import pprint
from opensbli.core.kernel import Kernel
from opensbli.core.bcs import BoundaryConditionBase


class split_bc(object):
    
    def __init__(self, *btypes):
        """Bc types are a list of tuples, with the first value in the tuple
        the order and the second value in the tuple actual boundary condition,
        this should be instantiated with palne==False"""
        self.boundary_condtions = {}
        for b in bctypes:
            if not isinstance(b[0], int):
                raise ValueError("The first argument in split bc should be a 
                    integer")
            if not isinstance(b[1], BoundaryConditionBase):
                raise ValueError("The second argument in split bc should be the
                    instantiated boundary condition")
            if b[1].full_plane:
                raise ValueError("The boundary condition %s should be instantiated
                                 with plane=False argument" % b[1])
            self.boundary_condtions[b[0]] = b[1]
        # The number of splits would be one less than the number of bc types
        self.nsplits = len(self.boundary_condtions.keys) - 1
        return

