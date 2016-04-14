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
    """ Class that defines data exchange on the same block. """
    def __init__(self, grid):
        range_of_evaluation = [tuple([0 + grid.halos[i][0], s + grid.halos[i][1]]) for i, s in enumerate(grid.shape)]
        # Size of transfers
        self.transfer_size = [r[1] - r[0] for r in range_of_evaluation]
        self.transfer_from = [grid.halos[i][0] for i, s in enumerate(grid.shape)]
        self.transfer_to = [grid.halos[i][0] for i, s in enumerate(grid.shape)]
        self.transfer_arrays = []
        return
class BoundaryClass(object):
    """
    Base class for boundary conditions, we store the name of the boundary condition and 
    type of the boundary for debugging purposes only.
    All the boundary conditions application requires this base class on the grid
    Computations can be computational Kernels or exchange type classes
    """
    def __init__(self, grid):
        self.boundaries = [None for sh in range(2) for sh in grid.shape]
        self.boundary_types = [None for sh in range(2) for sh in grid.shape]
        self.computations = [None for sh in range(2) for sh in grid.shape]
        return


class periodicboundary():
    """ class for applying periodic boundary condition, this updates the BoundaryClass specified"""
    def apply_boundary(self,boundaryclass, grid, arrays, boundary_direction, matching_face = None):
        if matching_face:
            raise NotImplementedError("Periodic boundary condition for different blocks")
        else:
            self.update_boundary_class(boundaryclass, boundary_direction, ['exchange_self', 'exchange_self'])
            left, right  = self.periodic_bc(boundary_direction, grid, arrays)
            boundaryclass.computations[boundary_direction*2 + 0] = left
            boundaryclass.computations[boundary_direction*2 + 1] = right
        return boundaryclass
    def update_boundary_class(self, boundaryclass, boundary_direction, btypes):
        # We donot require boundary type this is mainly for debugging
        boundaryclass.boundary_types[boundary_direction*2 + 0] = btypes[0]
        boundaryclass.boundary_types[boundary_direction*2 + 1] = btypes[1]

        return
    def periodic_bc(self, direction, grid, arrays, matching_face=None):
        """ Periodic boundary condition. """

        transfers = []
        # Generic transfer for the grid
        if matching_face:
            pass
        else:
            transfers_left = ExchangeSelf(grid)
            transfers_right = ExchangeSelf(grid)
        # Left transfers are from the start of grid to the end of grid (nx)
        transfers_left.transfer_size[direction] = abs(grid.halos[direction][0])
        transfers_left.transfer_from[direction] = 0
        transfers_left.transfer_to[direction] = grid.shape[direction]
        # Right transfers are from end of grid halo points to the start of the halo points
        transfers_right.transfer_size[direction] = abs(grid.halos[direction][0])
        transfers_right.transfer_from[direction] = grid.shape[direction]+ grid.halos[direction][0]
        transfers_right.transfer_to[direction] = grid.halos[direction][0]
        if matching_face:
            pass
        else:
            transfers_left.transfer_arrays = arrays
            transfers_right.transfer_arrays = arrays

        return transfers_left, transfers_right

class symmetry():
    """Applying symmetry boundary condition."""
    types = {0: 'Left', 1 :'Right'}
    def apply_boundary(self, boundaryclass, grid, arrays, boundary_direction, side):
        """
        input arrays should be a list of lists, vectors should be in the inner lists
        boundary_direction: the direction on the grid Symmetry Boundary condition should be applied
        side corresponds to the left or right face of boundary_direction
        """
        self.update_boundary_class(boundaryclass, boundary_direction, ['Computation', 'Computation'])
        boundaryclass.computations[boundary_direction*2 + side] = self.symmetry_bc(boundary_direction, side, grid, arrays)
        return boundaryclass
    
    def update_boundary_class(self, boundaryclass, boundary_direction, btypes):
        # We donot require boundary type this is mainly for debugging
        boundaryclass.boundary_types[boundary_direction*2 + 0] = btypes[0]
        boundaryclass.boundary_types[boundary_direction*2 + 1] = btypes[1]
        return
    
    def symmetry_bc(self, direction, side, grid, arrays):
        """ Writing symmetry bc as a kernel"""
        from .kernel import *
        if side == 0:
            base = 0 # Left side starting index
            
            tuples = [tuple([-t, t]) for t in range(1,abs(grid.halos[direction][side]) + 1)] # Indices to be updated relative to base
            
            # range of evaluation of the Kernel, Take the entire grid range and modify according to direction
            range_of_evaluation = [tuple([0 + grid.halos[direction][0], s + grid.halos[direction][1]]) for s in grid.shape]
            range_of_evaluation[direction] = tuple([base,base+1])
            
            # Get the equations for symmetry Bc
            symmetryeq = self.get_symmetry_equations(tuples, arrays,direction,grid)
            
            # Kernel for the computations
            kern = Kernel(symmetryeq, range_of_evaluation, "Symmetry bc %d %s"%(direction,self.types[side]) , grid)
            
            return kern
        
        elif side == 1:
            
            base = grid.shape[direction] # The end point of the domain in the direction of boundary
            
            tuples = [tuple([t,-t]) for t in range(1,abs(grid.halos[direction][side]) + 1)] # Indices to be updated relative to base
            
            # range of evaluation of the Kernel, First the entire grid range
            range_of_evaluation = [tuple([0 + grid.halos[direction][0], s + grid.halos[direction][1]]) for s in grid.shape]
            range_of_evaluation[direction] = tuple([base-1,base])
            
            # Get equations for the symmetry BC
            symmetryeq = self.get_symmetry_equations(tuples, arrays, direction,grid)
            
            # Kernel for the computations
            kern = Kernel(symmetryeq, range_of_evaluation, "Symmetry bc %d %s"%(direction,self.types[side]) , grid)
            
            return kern
        else:
            raise ValueError("The last input for symmetry boundary can be either\
                0 or 1 Corresponding to left or right boundary in the given direction")

        return
    def get_symmetry_equations(self, tuples, arrays, direction, grid):
        """ This returns the symmetry boundary condition equations depending on the direction and the 
        type of the variable
        Vector components in the 'direction' specified are reversed and the rest are kept the same
        Scalar is equated to the same value"""
        symmetryeq = []
        for ar in arrays:
            if isinstance(ar, list):
                """ Vector symmetry
                """
                newarray = [grid.work_array(a.base) for a in ar]
                newind = ar[0].indices
                direction_index = newind[direction]
                areq = []
                for tup in tuples:
                    lhs_arrays = [grid.make_array_grid(a.subs({direction_index:direction_index + tup[0]})) for\
                        number,a in enumerate(newarray)]
                    rhs_arrays = [grid.make_array_grid(a.subs({direction_index:direction_index + tup[1]})) \
                        if (number != direction) else \
                        -grid.make_array_grid(a.subs({direction_index:direction_index+tup[1]})) \
                            for number,a in enumerate(newarray)]
                    areq += [Eq(lh,rh, evaluate=False) for lh,rh in  zip(lhs_arrays,rhs_arrays)]
                symmetryeq += areq
            else:
                # Scalar symmetry
                areq = []
                newarray = grid.work_array(ar.base)
                direction_index = newarray.indices[direction]
                
                for tup in tuples:
                    lhs_arrays = [grid.make_array_grid(newarray.subs({direction_index:direction_index+tup[0]})) ]
                    rhs_arrays = [grid.make_array_grid(newarray.subs({direction_index:direction_index+tup[1]}))]
                    areq += [Eq(lh,rh, evaluate=False) for lh,rh in  zip(lhs_arrays,rhs_arrays)]
                
                symmetryeq += areq
                
        return symmetryeq
