from opensbli.core.boundary_conditions.bc_core import BoundaryConditionBase
from sympy import flatten, Matrix
from opensbli.core.kernel import Kernel
from opensbli.core.boundary_conditions.exchange import ExchangeSelf


class MultiBlockBoundary(object):
    pass


class InterfaceBC(BoundaryConditionBase, MultiBlockBoundary):
    def __init__(self, direction, side, match=(None, None, None), plane=True):
        # check if the match is a boundary type
        BoundaryConditionBase.__init__(self, direction, side, plane)
        self.match = match
        self.bc_name = "interface"
        return

    def apply(self, arrays, block):
        """"""
        return []

    def set_kernel_range(self, kernel, block):
        """ Sets the boundary condition kernel ranges based on direction and side.

        :arg object kernel: Computational boundary condition kernel.
        :arg object block: The SimulationBlock the boundary conditions are used on.
        :returns kernel: The computational kernel with updated ranges."""
        side, direction = self.side, self.direction
        kernel.ranges = block.ranges[:]
        # For interface boundary condition we should not include the boundary point
        if side == 0:
            left = 0
            right = 0
        elif side == 1:
            left = 0
            right = 0
        kernel.ranges[direction] = [block.ranges[direction][side]+left, block.ranges[direction][side]+right]
        return kernel

    def apply_interface(self, arrays, block, multiblock_descriptor, other_arrays=None):
        # halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        other_block = multiblock_descriptor.get_block(self.match[0])
        if other_arrays:
            other_block_arrays = other_arrays[:]
        else:
            other_block_arrays = [other_block.work_array(str(a.base.label)) for a in flatten(arrays)]

        # From corresponds to the block
        halos_block1 = self.get_halo_values(block)
        halos_block2 = self.get_halo_values(other_block)

        # Get the number of halos requred for block 2
        from_location = [d[0] for d in halos_block2]
        to_location = [d[0] for d in halos_block2]
        direction_block1, side_block1 = self.direction, self.side
        direction_block2, side_block2 = self.match[1], self.match[2]
        if self.full_plane is True:
            # No need to check
            idx_block1 = block.Idxed_shape
            idx_block2 = other_block.Idxed_shape
        # Set the start from depending on the side
        if side_block1 == 0:
            from_location[direction_block1] = idx_block1[direction_block1].lower
        else:
            from_location[direction_block1] = idx_block1[direction_block1].upper + halos_block2[direction_block2][0]
        # Check the side and modify transfer to of the other block
        if side_block2 == 0:
            # no need to modify the to starting point is the halos
            pass
        else:
            # to starting point is the end grid point
            to_location[direction_block2] = idx_block2[direction_block2].upper
        # size of transfers this would be dependant on block 2
        transfer_size = Matrix([i.upper + i.lower for i in idx_block2]) + \
            Matrix([abs(dire[0]) + abs(dire[1]) for dire in halos_block2])
        transfer_size[direction_block1] = abs(halos_block2[direction_block2][side_block2])
        ker = Kernel(block)
        # We will make use of exchange self currently, later we will change the perioidBC and combine them both
        exchange = ExchangeSelf(block, self.direction, self.side)
        exchange.computation_name = self.bc_name + "_exchange"
        exchange.set_transfer_size(transfer_size)
        exchange.set_transfer_from(from_location)
        exchange.set_transfer_to(to_location)
        # Manually set the arrays
        exchange.transfer_arrays = flatten(arrays)
        exchange.from_arrays = flatten(arrays)
        exchange.to_arrays = flatten(other_block_arrays)
        # ex.set_arrays(arrays)
        exchange.flip = self.match
        exchange.number = ker.kernel_no
        return exchange


class SharedInterfaceBC(InterfaceBC, BoundaryConditionBase, MultiBlockBoundary):

    def apply_interface(self, arrays, block, multiblock_descriptor, other_arrays=None):
        # halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        other_block = multiblock_descriptor.get_block(self.match[0])
        if other_arrays:
            other_block_arrays = other_arrays[:]
        else:
            other_block_arrays = [other_block.work_array(str(a.base.label)) for a in flatten(arrays)]

        # From corresponds to the block
        halos_block1 = self.get_halo_values(block)
        halos_block2 = self.get_halo_values(other_block)

        # Get the number of halos requred for block 2
        from_location = [d[0] for d in halos_block2]
        to_location = [d[0] for d in halos_block2]
        direction_block1, side_block1 = self.direction, self.side
        direction_block2, side_block2 = self.match[1], self.match[2]
        if self.full_plane is True:
            # No need to check
            idx_block1 = block.Idxed_shape
            idx_block2 = other_block.Idxed_shape
        # Set the start from depending on the side as the interface is shared we offset the point by 1 for min and -1 for max
        if side_block1 == 0:
            from_location[direction_block1] = idx_block1[direction_block1].lower + 1
        else:
            from_location[direction_block1] = idx_block1[direction_block1].upper + halos_block2[direction_block2][0] - 1
        # Check the side and modify transfer to of the other block
        if side_block2 == 0:
            # no need to modify the to starting point is the halos
            pass
        else:
            # to starting point is the end grid point
            to_location[direction_block2] = idx_block2[direction_block2].upper
        # size of transfers this would be dependant on block 2
        transfer_size = Matrix([i.upper + i.lower for i in idx_block2]) + \
            Matrix([abs(dire[0]) + abs(dire[1]) for dire in halos_block2])
        transfer_size[direction_block1] = abs(halos_block2[direction_block2][side_block2])
        ker = Kernel(block)
        # We will make use of exchange self currently, later we will change the perioidBC and combine them both
        exchange = ExchangeSelf(block, self.direction, self.side)
        exchange.computation_name = "Shared_interface" + "_exchange"
        exchange.set_transfer_size(transfer_size)
        exchange.set_transfer_from(from_location)
        exchange.set_transfer_to(to_location)
        # Manually set the arrays
        exchange.transfer_arrays = flatten(arrays)
        exchange.from_arrays = flatten(arrays)
        exchange.to_arrays = flatten(other_block_arrays)
        # ex.set_arrays(arrays)
        exchange.number = ker.kernel_no
        exchange.flip = self.match
        return exchange
