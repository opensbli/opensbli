from sympy import *

class Exchange(object):
    def __init__(self, grid):
        range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
        # Size of transfers
        self.transfer_size = [r[1] - r[0] for r in range_of_evaluation]
        self.transfer_from = [grid.halos[i][0] for i, s in enumerate(grid.shape)]
        self.transfer_to = [grid.halos[i][0] for i, s in enumerate(grid.shape)]
        self.transfer_arrays = []
        return


class BoundaryConditions(object):

    """ Boundary conditions applied to the equations by enforcing values at the boundary grid points. """

    types = {"periodic":"exchange_self", "symmetry":"exchange_self"}

    def __init__(self, bcs, grid, arrays):
        """ Initialise the boundary conditions. """

        if len(bcs) != len(grid.shape):
            raise ValueError("The number of boundary conditions and the dimensions of the grid do not match.")
        self.boundaries = bcs
        self.computations = [None for b in bcs for a in b]
        self.transfers = [None for b in bcs for a in b]
        self.type = self.get_type()
        for ind,bc in enumerate(self.boundaries):
            if bc[0] == bc[1] and bc[0] == "periodic":
                left, right = self.periodic_bc(ind, grid, arrays)
                self.transfers[ind*2 + 0] = left
                self.transfers[ind*2 + 1] = right
            else:
                raise NotImplementedError("Boundary condition %s not implemented" % bc)
        return

    def get_type(self):
        """ Return the type of boundary condition to be applied. """
        types = BoundaryConditions.types
        type_of_boundary = [[types[bc[0]], types[bc[1]]] for bc in self.boundaries]
        type_of_boundary = flatten(type_of_boundary)
        return type_of_boundary

    def periodic_bc(self, direction, grid, arrays):
        """ Periodic boundary condition. """

        transfers = []
        # Generic transfer for the grid
        transfers_left = Exchange(grid)
        transfers_right = Exchange(grid)
        # Left transfers are from the start of grid to the end of grid (nx)
        transfers_left.transfer_size[direction] = abs(grid.halos[direction][0])
        transfers_left.transfer_from[direction] = 0
        transfers_left.transfer_to[direction] = grid.shape[direction]
        transfers_left.transfer_arrays = arrays
        # Right transfers are from end of grid halo points to the start of the halo points
        transfers_right.transfer_size[direction] = abs(grid.halos[direction][0])
        transfers_right.transfer_from[direction] = grid.shape[direction]+ grid.halos[direction][0]
        transfers_right.transfer_to[direction] = grid.halos[direction][0]
        transfers_right.transfer_arrays = arrays
        return transfers_left, transfers_right


