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
        types = flatten(types)
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
        return


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

        kernel = Kernel(block, computation_name="Extrapolation boundary dir%d side%d" % (boundary_direction, side))
        kernel.set_boundary_plane_range(block, boundary_direction, side)
        halos = kernel.get_plane_halos(block)
        base = block.ranges[boundary_direction][side]
        if side == 0:
            from_side_factor = -1
            to_side_factor = 1
        elif side == 1:
            from_side_factor = 1
            to_side_factor = -1       
        indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[boundary_direction][side]) + 1)]
        equations = []
        for no, dset in enumerate(flatten(arrays)):
            array_equation = []
            loc = list(dset.indices)
            for index in indices:
                loc1, loc2 = loc[:], loc[:]
                loc1[boundary_direction] += index[0]
                loc2[boundary_direction] += index[1]
                array_equation += [Eq(dset.base[loc1], dset.base[loc2])]
            equations += array_equation
        kernel.add_equation(equations)
        kernel.set_boundary_plane_range(block, boundary_direction, side) 
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
        kernel.set_halo_range(boundary_direction, side, block.boundary_halos[boundary_direction][side])
        kernel.add_equation(self.equations)
        kernel.update_block_datasets(block)
        return kernel

class SymmetryBoundaryConditionBlock(BoundaryConditionBase):

    def __init__(self, boundary_direction, side, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        return

    def apply(self, arrays, boundary_direction, side, block):
        metric = zeros(block.ndim)
        for i in range(block.ndim):
            metric[i,i] = DataObject('D%d%d' % (i,i))
        for row in range(block.ndim):
            total = []
            expr = metric[row, :]
            for item in expr:
                total += [item**2]
            metric[row, :] = metric[row, :] / sqrt(sum(total))

        lhs_eqns = flatten(arrays)
        rhs_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                transformed_vector = metric*Matrix(ar)
                transformed_vector[boundary_direction] = -1*transformed_vector[boundary_direction]
                rhs_eqns += flatten(transformed_vector)
            else:
                rhs_eqns += [ar]

        kernel = Kernel(block, computation_name="Symmetry boundary dir%d side%d" % (boundary_direction, side))
        kernel.set_boundary_plane_range(block, boundary_direction, side)
        halos = kernel.get_plane_halos(block)      
        base = block.ranges[boundary_direction][side]
       
        if side == 0:
            from_side_factor = -1
            to_side_factor = 1
        elif side == 1:
            from_side_factor = 1
            to_side_factor = -1
                    
        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[boundary_direction][side]) + 1)]
        loc = list(lhs_eqns[0].indices)
        final_equations = []
        for index in transfer_indices:
            array_equations = []
            loc_lhs, loc_rhs = loc[:], loc[:]
            loc_lhs[boundary_direction] += index[0]
            loc_rhs[boundary_direction] += index[1]
            for left, right in zip(lhs_eqns, rhs_eqns):
                left = self.convert_dataset_base_expr_to_datasets(left, loc_lhs)
                right = self.convert_dataset_base_expr_to_datasets(right, loc_rhs)
                array_equations += [Eq(left, right, evaluate=False)]
            pprint(array_equations)
            final_equations += array_equations
        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)

        return kernel

    def convert_dataset_base_expr_to_datasets(self, expression, index):
        for a in expression.atoms(DataSet):
            b = a.base
            expression = expression.xreplace({a: b[index]})
        return expression



