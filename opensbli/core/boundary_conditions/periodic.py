from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase
from opensbli.core.boundary_conditions.exchange import ExchangeSelf
from opensbli.core.kernel import Kernel
from sympy import Matrix


class PeriodicBC(BoundaryConditionBase):
    """ Applies an exchange periodic boundary condition.

    :arg int direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, direction, side, plane=True):
        BoundaryConditionBase.__init__(self, direction, side, plane)
        return

    def halos(self):
        return True

    def apply(self, arrays, block):
        # Get the exchanges which form the computations.
        if self.full_plane:
            exchange = self.get_exchange_plane(arrays, block)
        return exchange

    def get_exchange_plane(self, arrays, block):
        """ Create the exchange computations which copy the block point values to/from the periodic domain boundaries. """

        # Create a kernel this is a neater way to implement the transfers
        ker = Kernel(block)
        halos = self.get_halo_values(block)
        size, from_location, to_location = self.get_transfers(block.Idxed_shape, halos)
        ex = ExchangeSelf(block, self.direction, self.side)
        ex.set_transfer_size(size)
        ex.set_transfer_from(from_location)
        ex.set_transfer_to(to_location)
        ex.set_arrays(arrays)
        ex.number = ker.kernel_no
        return ex

    def get_transfers(self, idx, halos):
        direction, side = self.direction, self.side
        transfer_from = [d[0] for d in halos]
        transfer_to = [d[0] for d in halos]
        if side == 0:
            transfer_from[direction] = idx[direction].lower
            transfer_to[direction] = idx[direction].upper
        else:
            transfer_from[direction] = idx[direction].upper + halos[direction][0]
            transfer_to[direction] = idx[direction].lower + halos[direction][0]

        transfer_size = Matrix([i.upper + i.lower for i in idx]) + \
            Matrix([abs(dire[0]) + abs(dire[1]) for dire in halos])
        transfer_size[direction] = abs(halos[direction][side])
        return transfer_size, transfer_from, transfer_to
