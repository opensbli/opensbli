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

    def check_modify_central(self):
        modify = {}
        for no, val in enumerate(self.boundary_types):
            left = val[0]; right = val[1]
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

from sympy.functions.elementary.piecewise import ExprCondPair
class Carpenter(object):
    def __init__(self):

        self.bc4_coefficients = self.carp4_coefficients()
        self.bc4_symbols = self.carp4_symbols()
        self.bc4_2_symbols = self.second_der_symbols()
        self.bc4_2_coefficients = self.second_der_coefficients()
        return

    def function_points(self, expression, direction, side):
        f_matrix = zeros(6,6)
        loc = list(list(expression.atoms(DataSet))[0].indices)
        for shift in range(6):
            func_points = []
            for index in range(6):
                new_loc = loc[:]
                new_loc[direction] += index - shift
                for dset in expression.atoms(DataSet):
                    expression = expression.replace(dset, dset.base[new_loc])
                func_points.append(expression)
            f_matrix[:,shift] = func_points
        if side == 0:
            f_matrix = f_matrix[:,0:4]
        elif side == 1:
            print "this should be negative or modified according to sbli, check"
            f_matrix = f_matrix.transpose()[:,0:4]
        else:
            raise NotImplementedError("Side must be 0 or 1")
        return f_matrix

    def weight_function_points(self, func_points, direction, order, block, char_BC=False):
        # WARNING side == 1, take negative h value? Check SBLI

        if order == 1:
            h = (block.deltas[direction])**(-1)
            if char_BC:
                weighted = h*(self.bc4_symbols[0,:]*func_points)
            else:
                weighted = zeros(4,1)
                for i in range(4):
                    weighted[i] = h*(self.bc4_coefficients[i,:]*func_points[:,i])
        elif order == 2:
            h_sq = (block.deltas[direction])**(-2)
            weighted = zeros(2,1)
            for i in range(2):
                weighted[i] = h_sq*(self.bc4_2_coefficients[i,:]*func_points[0:5,i])
        else:
            raise NotImplementedError("Only 1st and 2nd derivatives implemented")
        return weighted

    def expr_cond_pairs(self, fn, direction, side, order, block):
        fn_pts = self.function_points(fn, direction, side)
        derivatives = self.weight_function_points(fn_pts, direction, order, block)
        idx = block.grid_indexes[direction]
        if side == 0:
            mul_factor = 1
            start = block.ranges[direction][side]
        else:
            mul_factor = -1
            start = block.ranges[direction][side] - 1
        ecs = []
        for no,d in enumerate(derivatives):
            loc = start + mul_factor*no
            ecs += [ExprCondPair(d, Equality(idx,loc))]
        return ecs

    def carp4_symbols(self):
        """
        Function to return symbolic bc4 matrix for testing
        """
        bc4 = MatrixSymbol('BC', 4, 6)
        bc4 = Matrix(bc4)
        return bc4

    def second_der_symbols(self):
        bc4_2 = MatrixSymbol('BCC', 2, 5)
        bc4_2 = Matrix(bc4_2)
        return bc4_2

    def second_der_coefficients(self):
        bc4_2 = Matrix([[35.0, -104.0, 114.0, -56.0, 11.0], [11.0, -20.0, 6.0, 4.0, -1.0]])/12.0
        for i in range(bc4_2.shape[0]):
            for j in range(bc4_2.shape[1]):
                bc4_2[i,j] = nsimplify(bc4_2[i,j])
        return bc4_2

    def carp4_coefficients(self):
        """
        Function to return the bc4 matrix containing coefficients for the Carpenter wall boundary treatment.
        """
        R1 = -(2177.0*sqrt(295369.0)-1166427.0)/25488.0
        R2 = (66195.0*sqrt(53.0)*sqrt(5573.0)-35909375.0)/101952.0

        al4_0 = [-(216.0*R2+2160.0*R1-2125.0)/12960.0, (81.0*R2+675.0*R1+415.0)/540.0, -(72.0*R2+720.0*R1+445.0)/1440.0, -(108.0*R2+756.0*R1+421.0)/1296.0]
        al4_1 = [(81.0*R2+675.0*R1+415.0)/540.0, -(4104.0*R2+32400.0*R1+11225.0)/4320.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(216.0*R2+2160.0*R1+655.0)/4320.0]
        al4_2 = [-(72.0*R2+720.0*R1+445.0)/1440.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(4104.0*R2+32400.0*R1+12785.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0]
        al4_3 = [-(108.0*R2+756.0*R1+421.0)/1296.0, -(216.0*R2+2160.0*R1+655.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0, -(216.0*R2+2160.0*R1-12085.0)/12960.0]

        al4 = Matrix([al4_0,al4_1,al4_2,al4_3])

        ar4_0 = [(-1.0)/2.0, -(864.0*R2+6480.0*R1+305.0)/4320.0, (216.0*R2+1620.0*R1+725.0)/540.0,-(864.0*R2+6480.0*R1+3335.0)/4320.0, 0.0, 0.0]
        ar4_1 = [(864.0*R2+6480.0*R1+305.0)/4320.0, 0.0, -(864.0*R2+6480.0*R1+2315.0)/1440.0, (108.0*R2+810.0*R1+415.0)/270.0, 0.0, 0.0]
        ar4_2 = [-(216.0*R2+1620.0*R1+725.0)/540.0, (864.0*R2+6480.0*R1+2315.0)/1440.0, 0.0, -(864.0*R2+6480.0*R1+785.0)/4320.0, -1.0/12.0, 0.0]
        ar4_3 = [(864.0*R2+6480.0*R1+3335.0)/4320.0, -(108.0*R2+810.0*R1+415.0)/270.0, (864.0*R2+6480.0*R1+785.0)/4320.0, 0.0, 8.0/12.0, -1.0/12.0]
        ar4 = Matrix([ar4_0, ar4_1, ar4_2, ar4_3])

        # Form inverse and convert to rational
        al4_inv = al4.inv()
        bc4 = al4_inv*ar4
        for i in range(bc4.shape[0]):
            for j in range(bc4.shape[1]):
                bc4[i,j] = nsimplify(bc4[i,j])
        return bc4

class ModifyCentralDerivative(object):
    """ A place holder for the boundary conditions on which the central derivative should be modified"""
    pass

class AdiabaticWall(ModifyCentralDerivative, BoundaryConditionBase):
    def __init__(self, boundary_direction, side, scheme = None,plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return
    def apply(self,  arrays, boundary_direction, side, block):

        print "Implement the scheme for adiabatic wall boundary condition"
        # Testing the boundary condition with a linear extrapolation
        kernel = LinearExtrapolateBoundaryConditionBlock(boundary_direction, side, plane=True).apply(arrays, boundary_direction, side, block)
        return kernel

def apply_modify_derivative(order, fn, bcs, block, value):

    return
