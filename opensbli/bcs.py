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

from sympy import *


class ExchangeSelf(object):

    """ Defines data exchange on the same block. """

    def __init__(self, grid):
        # Range of evaluation (i.e. the grid points, including the halo points, over which the computation should be performed).
        range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
        # Size of transfers
        self.transfer_size = [r[1] - r[0] for r in range_of_evaluation]
        self.transfer_from = [grid.halos[i][0] for i, s in enumerate(grid.shape)]
        self.transfer_to = [grid.halos[i][0] for i, s in enumerate(grid.shape)]
        self.transfer_arrays = []
        return


class BoundaryConditionBase(object):

    """ Base class for boundary conditions. We store the name of the boundary condition and type of the boundary for debugging purposes only.
    The application of the boundary conditions requires this base class on the grid.
    Computations can be computational Kernels or Exchange type objects.
    """

    def __init__(self, grid):
        self.grid = grid
        self.boundary_types = [None for sh in range(2) for sh in grid.shape]
        self.computations = [None for sh in range(2) for sh in grid.shape]
        return

class PeriodicBoundaryCondition(BoundaryConditionBase):

    """ Periodic boundary condition. This updates the BoundaryClass specified"""

    def apply(self, arrays, boundary_direction, matching_face=None):
        if matching_face:
            raise NotImplementedError("Periodic boundary condition for different blocks.")
        else:
            # Set boundary type. We do not really require the boundary type; this is mainly for debugging.
            self.boundary_types[boundary_direction*2 + 0] = 'exchange_self'
            self.boundary_types[boundary_direction*2 + 1] = 'exchange_self'
            
            # Get the exchanges which form the computations.
            left, right = self.get_exchange(boundary_direction, arrays)
            self.computations[boundary_direction*2 + 0] = left
            self.computations[boundary_direction*2 + 1] = right
        return

    def get_exchange(self, direction, arrays, matching_face=None):
        """ Create the exchange computations which copy the grid point values to/from the periodic domain boundaries. """

        # Generic transfer for the grid
        if matching_face:
            pass
        else:
            transfers_left = ExchangeSelf(self.grid)
            transfers_right = ExchangeSelf(self.grid)

        # Left transfers are from the start of grid to the end of grid (nx)
        transfers_left.transfer_size[direction] = abs(self.grid.halos[direction][0])
        transfers_left.transfer_from[direction] = 0
        transfers_left.transfer_to[direction] = self.grid.shape[direction]

        # Right transfers are from end of grid halo points to the start of the halo points
        transfers_right.transfer_size[direction] = abs(self.grid.halos[direction][0])
        transfers_right.transfer_from[direction] = self.grid.shape[direction] + self.grid.halos[direction][0]
        transfers_right.transfer_to[direction] = self.grid.halos[direction][0]
        if matching_face:
            pass
        else:
            transfers_left.transfer_arrays = arrays
            transfers_right.transfer_arrays = arrays

        return transfers_left, transfers_right


class SymmetryBoundaryCondition(BoundaryConditionBase):

    """ Symmetry boundary condition. """
    types = {0: 'Left', 1: 'Right'}

    def apply(self, arrays, boundary_direction, side):
        """ Apply the symmetry boundary condition.
        
        :arg grid: The grid on which the boundary condition is applied.
        :arg arrays: A list of lists. vectors should be in the inner lists.
        :arg boundary_direction: The direction on the grid symmetry boundary condition should be applied.
        :arg side: Corresponds to the left or right face of boundary_direction.
        """
        # Set boundary type. We do not really require the boundary type; this is mainly for debugging.
        self.boundary_types[boundary_direction*2 + 0] = 'Computation'
        self.boundary_types[boundary_direction*2 + 1] = 'Computation'
        
        # Create the kernel.
        self.computations[boundary_direction*2 + side] = self.get_kernel(boundary_direction, side, arrays)
        return

    def get_kernel(self, direction, side, arrays):
        """ Write the application of the symmetry boundary condition as a Kernel. """
        from .kernel import *
        if side == 0:
            base = 0  # Left side starting index

            tuples = [tuple([-t, t]) for t in range(1, abs(self.grid.halos[direction][side]) + 1)]  # Indices to be updated relative to base

            # range of evaluation of the Kernel, Take the entire grid range and modify according to direction
            range_of_evaluation = [tuple([0 + self.grid.halos[direction][0], s + self.grid.halos[direction][1]]) for s in self.grid.shape]
            range_of_evaluation[direction] = tuple([base, base+1])

            # Get the equations for symmetry Bc
            symmetry_equations = self.get_symmetry_equations(tuples, arrays, direction)

            # Kernel for the computations
            kernel = Kernel(symmetry_equations, range_of_evaluation, "Symmetry bc %d %s" % (direction, self.types[side]), self.grid)

            return kernel

        elif side == 1:

            base = self.grid.shape[direction]  # The end point of the domain in the direction of the boundary

            tuples = [tuple([t, -t]) for t in range(1, abs(self.grid.halos[direction][side]) + 1)]  # Indices to be updated relative to base

            # range of evaluation of the Kernel, First the entire grid range
            range_of_evaluation = [tuple([0 + self.grid.halos[direction][0], s + self.grid.halos[direction][1]]) for s in self.grid.shape]
            range_of_evaluation[direction] = tuple([base-1, base])

            # Get equations for the symmetry BC
            symmetry_equations = self.get_symmetry_equations(tuples, arrays, direction)

            # Kernel for the computations
            kernel = Kernel(symmetry_equations, range_of_evaluation, "Symmetry bc %d %s" % (direction, self.types[side]), self.grid)

            return kernel
            
        else:
            raise ValueError("The 'side' of the symmetry boundary should be either 0 or 1, corresponding to left or right boundary in the given direction.")


    def get_symmetry_equations(self, tuples, arrays, direction):
        """ Return the symmetry boundary condition equations depending on the direction and the type of the variable.
        Vector components in the 'direction' specified are reversed and the rest are kept the same.
        Scalar is equated to the same value. """
        symmetry_equations = []
        for array in arrays:
            if isinstance(array, list):
                # Vector symmetry
                new_array = [self.grid.work_array(a.base) for a in array]
                new_indices = array[0].indices
                direction_index = new_indices[direction]
                array_equation = []
                for t in tuples:
                    lhs_arrays = []
                    for number, a in enumerate(new_array):
                        lhs_arrays.append(self.grid.get_array_on_grid(a.subs({direction_index: direction_index + t[0]})))

                    rhs_arrays = []
                    for number, a in enumerate(new_array):
                        if (number != direction):
                            rhs_arrays.append(self.grid.get_array_on_grid(a.subs({direction_index: direction_index + t[1]})))
                        else:
                            rhs_arrays.append(-self.grid.get_array_on_grid(a.subs({direction_index: direction_index + t[1]})))

                    array_equation += [Eq(lhs, rhs, evaluate=False) for lhs, rhs in zip(lhs_arrays, rhs_arrays)]
                symmetry_equations += array_equation
            else:
                # Scalar symmetry
                array_equation = []
                new_array = self.grid.work_array(array.base)
                direction_index = new_array.indices[direction]

                for t in tuples:
                    lhs_arrays = [self.grid.get_array_on_grid(new_array.subs({direction_index: direction_index + t[0]}))]
                    rhs_arrays = [self.grid.get_array_on_grid(new_array.subs({direction_index: direction_index + t[1]}))]
                    array_equation += [Eq(lhs, rhs, evaluate=False) for lhs, rhs in zip(lhs_arrays, rhs_arrays)]

                symmetry_equations += array_equation

        return symmetry_equations
