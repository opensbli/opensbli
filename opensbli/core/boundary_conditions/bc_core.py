"""@brief Implementation of boundary conditions available in OpenSBLI.
   @authors David J Lusher, Satya Pramod Jammy
   @contributors
   @details
"""
from sympy import flatten, Idx
from opensbli.core.kernel import Kernel, ConstantsToDeclare
from opensbli.core.opensbliobjects import DataSet, ConstantIndexed
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.datatypes import Int
from opensbli.utilities.helperfunctions import get_min_max_halo_values
side_names = {0: 'left', 1: 'right'}


class ModifyCentralDerivative(object):
    """ A place holder for the boundary conditions on which the central derivative should be modified"""
    pass


class BoundaryConditionTypes(object):
    """ Base class for boundary condition types. We store the name of the boundary condition and type of the boundary for debugging purposes only.
    The application of the boundary conditions requires this base class on the grid.
    Computations can be computational Kernels or Exchange type objects."""

    def set_boundary_types(self, types, block):
        """ Adds the boundary types of the grid """
        # convert the list of types into a list of tuples
        types = flatten(types)
        self.check_boundarysizes_ndim_match(types)
        for t in types:
            t.convert_dataobject_to_dataset(block)
        it = iter(types)
        self.boundary_types = list(zip(it, it))
        return

    def check_boundarysizes_ndim_match(self, types):
        if len(types) != self.ndim*2:
            raise ValueError("Boundaries provided should match the number of dimension")
        return

    def check_modify_central(self):
        modify = {}
        for no, val in enumerate(self.boundary_types):
            left = val[0]
            right = val[1]
            if isinstance(left, ModifyCentralDerivative):
                if no in modify:
                    modify[no][0] = left
                else:
                    modify[no] = [left, None]
            if isinstance(right, ModifyCentralDerivative):
                if no in modify:
                    modify[no][1] = right
                else:
                    modify[no] = [None, right]
        return modify


class BoundaryConditionBase(object):
    """ Base class for common functionality between all boundary conditions.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, plane):
        if plane:
            self.full_plane = True
        else:
            self.full_plane = False
        self.direction = boundary_direction
        self.side = side
        self.equations = None
        return

    def convert_dataobject_to_dataset(self, block):
        """ Converts DataObjects to DataSets.
        :arg object block: OpenSBLI SimulationBlock."""
        if isinstance(self, SplitBC):
            for bc in self.bc_types:
                if bc.equations:
                    bc.equations = block.dataobjects_to_datasets_on_block(bc.equations)
        else:
            if self.equations:
                self.equations = block.dataobjects_to_datasets_on_block(self.equations)
        return

    def convert_dataset_base_expr_to_datasets(self, expression, index):
        """ Converts an expression containing DataSetBases to Datasets and updates locations.

        :arg object expression: Symbolic expression.
        :arg int index: Index to increment the DataSet by.
        :returns: object: expression: Updated symbolic expression."""
        for a in expression.atoms(DataSet):
            b = a.base
            expression = expression.xreplace({a: b[index]})
        return expression

    def generate_boundary_kernel(self, block, bc_name):
        if self.full_plane:
            return self.bc_plane_kernel(block, bc_name)
        else:
            return self.arbitrary_bc_plane_kernel(block, bc_name)

    def arbitrary_bc_plane_kernel(self, block, bc_name):
        """ Creates a computational kernel for use with Split BC."""
        bc_name = self.bc_name
        direction, side, split_number = self.direction, self.side, self.split_number
        kernel = Kernel(block, computation_name="%s bc direction-%d side-%d split-%d" % (bc_name, direction, side, split_number))
        print(kernel.computation_name)
        numbers = Idx('no', 2*block.ndim)
        ranges = ConstantIndexed('split_range_%d%d%d' % (direction, side, split_number), numbers)
        ranges.datatype = Int()
        kernel.ranges = ranges
        halo_ranges = ConstantIndexed('split_halo_range_%d%d%d' % (direction, side, split_number), numbers)
        halo_ranges.datatype = Int()
        kernel.halo_ranges = halo_ranges
        ConstantsToDeclare.add_constant(ranges)
        ConstantsToDeclare.add_constant(halo_ranges)
        halo_values = self.get_halo_values(block)
        return halo_values, kernel

    def set_kernel_range(self, kernel, block):
        """ Sets the boundary condition kernel ranges based on direction and side.

        :arg object kernel: Computational boundary condition kernel.
        :arg object block: The SimulationBlock the boundary conditions are used on.
        :returns kernel: The computational kernel with updated ranges."""
        side, direction = self.side, self.direction
        kernel.ranges = block.ranges[:]
        if side == 0:
            left = 0
            right = 1
        elif side == 1:
            left = -1
            right = 0
        kernel.ranges[direction] = [block.ranges[direction][side]+left, block.ranges[direction][side]+right]
        return kernel

    def get_halo_values(self, block):
        """ Gets the maximum numerical halo values.

        :arg object block: The SimulationBlock the boundary conditions are used on.
        :returns halo_values: Numerical values of the halos in all directions."""
        halo_values = []
        halo_objects = block.boundary_halos
        for i in range(len(halo_objects)):
            halo_m, halo_p = get_min_max_halo_values(halo_objects)
            halo_m, halo_p = halo_m[0], halo_p[0]
            halo_values.append([halo_m, halo_p])
        return halo_values

    def bc_plane_kernel(self, block, bc_name):
        direction, side = self.direction, self.side
        kernel = Kernel(block, computation_name="%s boundary dir%d side%d" % (bc_name, direction, side))
        kernel = self.set_kernel_range(kernel, block)
        halo_values = self.get_halo_values(block)
        # Add the halos to the kernel in directions not equal to boundary direction
        for i in [x for x in range(block.ndim) if x != direction]:
            kernel.halo_ranges[i][0] = block.boundary_halos[i][0]
            kernel.halo_ranges[i][1] = block.boundary_halos[i][1]
        return halo_values, kernel

    def create_boundary_equations(self, left_arrays, right_arrays, transfer_indices):
        """ Creates boundary equations for the given indices."""
        direction = self.direction
        if isinstance(left_arrays, list):
            loc = list(left_arrays[0].indices)
        else:
            loc = left_arrays.indices
        final_equations = []
        for index in transfer_indices:
            array_equations = []
            loc_lhs, loc_rhs = loc[:], loc[:]
            loc_lhs[direction] += index[0]
            loc_rhs[direction] += index[1]
            for left, right in zip(left_arrays, right_arrays):
                left = self.convert_dataset_base_expr_to_datasets(left, loc_lhs)
                right = self.convert_dataset_base_expr_to_datasets(right, loc_rhs)
                array_equations += [OpenSBLIEq(left, right, evaluate=False)]
            final_equations += array_equations
        return final_equations

    def set_side_factor(self):
        """ Sets the +/- 1 side factors for boundary condition halo numbering."""
        if self.side == 0:
            from_side_factor = -1
            to_side_factor = 1
        elif self.side == 1:
            from_side_factor = 1
            to_side_factor = -1
        return from_side_factor, to_side_factor


class SplitBC(BoundaryConditionBase):
    """ Functionality to apply more then one boundary condition along a boundary."""

    def __init__(self, boundary_direction, side, bcs, plane=False):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_types = bcs
        return

    # def pre_process_bc(self):
    #     return

    def apply(self, arrays, block):
        kernels = []
        for no, bc in enumerate(self.bc_types):
            bc.full_plane = False
            bc.split_number = no
            kernels.append(bc.apply(arrays, block))
        return kernels
