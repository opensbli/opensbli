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

# @author: New structure implemented by Satya P Jammy (October, 2016)

from sympy import *
from .kernel import *
side_names  = {0:'left', 1:'right'}
class Exchange(object):
    pass
class ExchangeSelf(Exchange):

    """ Defines data exchange on the same block. """

    def __init__(self, block, direction, side):
        ## Range of evaluation (i.e. the grid points, including the halo points, over which the computation should be performed).
        self.computation_name = "exchange"
        self.block_number = block.blocknumber
        self.block_name = block.blockname
        self.direction = direction
        self.side = side_names[side]
        return
    @property
    def name(self):
        return "%s%d"%(self.computation_name, self.number)

    def set_arrays(self, arrays):
        self.transfer_arrays = flatten(arrays)
        self.from_arrays = flatten(arrays)
        self.to_arrays = flatten(arrays)
        return

    def set_transfer_from(self, transfer):
        self.transfer_from = transfer
        return

    def set_transfer_to(self, transfer):
        self.transfer_to = transfer
        return
    def set_transfer_size(self, size):
        self.transfer_size = size
        return
    @property
    def algorithm_node_name(self):
        name = "Boundary_exchange_block_%d_direction%d_side%s"%(self.block_number,self.direction, self.side)
        return name

    def write_latex(self, latex):
        latex.write_string("This is an exchange self kernel on variables %s\\\\"%', '.join([str(a) for a in self.transfer_arrays]))
        return
    @property
    def opsc_code(self):
        """The string for calling the boundary condition in OPSC is update while creating 
        the code for exchanging data
        """
        return [self.call_name]

    #def set

class Wall(object):
    """Main class for wall boundary condition all wall bc types should be derived from here"""
    def halos(self):
        return False

class BoundaryConditionTypes(object):

    """ Base class for boundary conditions. We store the name of the boundary condition and type of the boundary for debugging purposes only.
    The application of the boundary conditions requires this base class on the grid.
    Computations can be computational Kernels or Exchange type objects.
    """

    def set_boundary_types(self, types):
        """ Adds the boundary types of the grid """
        # convert the list of types into a list of tuples
        self.check_boundarysizes_ndim_match(types)
        it = iter(types)
        self.boundary_types = zip(it, it)
        return

    def check_boundarysizes_ndim_match(self, types):
        if len(types) != self.ndim*2:
            raise ValueError("Boundaries provided should match the number of dimension")
        return

    def checkwallboundary(self):
        wallbcs = []
        for no, val in self.boundary_types:
            if isinstance(val, Wall):
                wallbcs += [no]
        return
class BoundaryConditionBase(object):

    def __init__(self, boundary_direction, side, plane):
        if plane:
            self.full_plane = True
        else:
            raise ValueError("")
        self.direction = boundary_direction
        self.side = side
        self.equations = None


class PeriodicBoundaryConditionBlock(BoundaryConditionBase):
    """
    Applies periodic boundary condition
    """
    def __init__(self, boundary_direction, side, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        return
    def halos(self):
        return True
    def apply(self, arrays, boundary_direction, side, block):
        # Get the exchanges which form the computations.
        exchange = self.get_exchange(arrays, boundary_direction,side, block)
        return exchange

    def get_exchange(self, arrays, boundary_direction, side, block):
        """ Create the exchange computations which copy the block point values to/from the periodic domain boundaries. """

        # Create a kernel this is a neater way to implement the transfers
        ker = Kernel(block)
        halos = ker.get_plane_halos(block)
        transfer_from = [halos[d][0] for d in range(block.ndim)]
        transfer_to = [halos[d][0] for d in range(block.ndim)]
        if side == 0:
            transfer_from[boundary_direction] = block.Idxed_shape[boundary_direction].lower
            transfer_to[boundary_direction] = block.Idxed_shape[boundary_direction].upper
        else:
            transfer_from[boundary_direction] = block.Idxed_shape[boundary_direction].upper + halos[boundary_direction][0]
            transfer_to[boundary_direction] = block.Idxed_shape[boundary_direction].lower + halos[boundary_direction][0]
        transfer_size = [block.Idxed_shape[dire].upper + block.Idxed_shape[dire].lower + abs(halos[dire][0]) + abs(halos[dire][1]) for dire in range(block.ndim)]
        transfer_size[boundary_direction] = abs(halos[boundary_direction][side])
        ex = ExchangeSelf(block, boundary_direction, side)
        ex.set_transfer_size(transfer_size)
        ex.set_transfer_from(transfer_from)
        ex.set_transfer_to(transfer_to)
        ex.set_arrays(arrays)
        ex.number = ker.kernel_no
        #exit()
        return ex

class SymmetryBoundaryCondition(object):
    pass

class LinearExtrapolateBoundaryConditionBlock(BoundaryConditionBase):
    """
    Applies zero gradient linear extrapolation boundary condition.
    """
    def __init__(self, boundary_direction, side, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        return
    def halos(self):
        return True
    def apply(self, arrays, boundary_direction, side, block):
        dire = boundary_direction
        kernel = Kernel(block, computation_name="Extrapolate boundary dir%d side%d" % (dire, side))
        halos = kernel.get_plane_halos(block)

        if side == 0:
            base = 0
            indices = [tuple([-t, t]) for t in range(1, abs(halos[dire][side]) + 1)]
            range_of_evaluation = [tuple([0 + halos[dire][0], s + halos[dire][1]]) for s in block.shape]
            range_of_evaluation[dire] = tuple([base, base+1])

        elif side == 1:
            base = block.shape[dire]  # The end point of the domain in the direction of the boundary
            indices = [tuple([t, -t]) for t in range(1, abs(halos[dire][side]) + 1)]
            range_of_evaluation = [tuple([0 + halos[dire][0], s + halos[dire][1]]) for s in block.shape]
            range_of_evaluation[dire] = tuple([base-1, base])
        # print "indices are: ", indices
        # print "range of eval is: ", range_of_evaluation

        equations = []
        for no, dset in enumerate(flatten(arrays)):
            array_equation = []
            loc = list(dset.indices)
            for idx in indices:
                loc1, loc2 = loc[:], loc[:]
                loc1[dire] += idx[0]
                loc2[dire] += idx[1]
                array_equation += [Eq(dset.base[loc1], dset.base[loc2])]
            equations += array_equation
        kernel.add_equation(equations)
        kernel.set_boundary_plane_range(block, dire, side) 
        kernel.update_block_datasets(block)
        return kernel


class DirichletBoundaryConditionBlock(BoundaryConditionBase):
    """Applies constant value Dirichlet boundary condition."""
    def __init__(self, boundary_direction, side , equations, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.equations = equations
        return
    def halos(self):
        return True
    def apply(self, arrays, boundary_direction, side, block):
        kernel = Kernel(block, computation_name="Dirichlet boundary dir%d side%d" % (boundary_direction, side))
        kernel.set_boundary_plane_range(block, boundary_direction, side)
        kernel.ranges[boundary_direction][side] = block.ranges[boundary_direction][side]
        # kernel.set_boundary_plane_range(block, boundary_direction, side) 
        print kernel.ranges, "string"
        # exit()
        kernel.set_halo_range(boundary_direction, side, block.boundary_halos[boundary_direction][side])
        # print kernel.ranges
        # print "halos are: ", halos
        # exit()
        # if side == 0:
        #     base = 0  # Left side starting index
        #     range_of_evaluation = [tuple([0 + halos[dire][0], s + halos[dire][1]]) for s in block.shape]
        #     range_of_evaluation[dire] = tuple([halos[dire][0], base])
        #     print "side is: ", side
        #     print "range is: ", range_of_evaluation
        #     print "\n"
        # elif side == 1:
        #     base = block.shape[dire]
        #     range_of_evaluation = [tuple([0 + halos[dire][0], s + halos[dire][1]]) for s in block.shape]
        #     range_of_evaluation[dire] = tuple([base, base + halos[dire][1]])
        #     print "side is: ", side
        #     print "range is: ", range_of_evaluation
        #     print "\n"
        # # eqns = [Eq(x,y) for x,y in zip(flatten(arrays), self.equations)]
        kernel.add_equation(self.equations)
        # kernel.set_boundary_plane_range(block, dire, side) 
        kernel.update_block_datasets(block)
        return kernel




