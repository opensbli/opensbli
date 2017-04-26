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
from opensbli.physical_models.euler_eigensystem import EulerEquations
from sympy import *
from .kernel import *
from opensbli.utilities.helperfunctions import increment_dataset
side_names  = {0:'left', 1:'right'}

class Exchange(object):
    pass

class ModifyCentralDerivative(object):
    """ A place holder for the boundary conditions on which the central derivative should be modified"""
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

    def convert_dataset_base_expr_to_datasets(self, expression, index):
        for a in expression.atoms(DataSet):
            b = a.base
            expression = expression.xreplace({a: b[index]})
        return expression

    def generate_boundary_kernel(self, direction, side, block, bc_name):
        kernel = Kernel(block, computation_name="%s boundary dir%d side%d" % (bc_name, direction, side))
        kernel.set_boundary_plane_range(block, direction, side)
        halos = kernel.get_plane_halos(block)
        # Add the halos to the kernel in directions not equal to boundary direction
        for i in [x for x in range(block.ndim) if x != direction]:
            kernel.halo_ranges[i][0] = block.boundary_halos[i][0]
            kernel.halo_ranges[i][1] = block.boundary_halos[i][1]
        return halos, kernel

    def create_boundary_equations(self, direction, left_arrays, right_arrays, transfer_indices):
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
                array_equations += [Eq(left, right, evaluate=False)]
            final_equations += array_equations
        return final_equations

    def set_side_factor(self, side):
        if side == 0:
            from_side_factor = -1
            to_side_factor = 1
        elif side == 1:
            from_side_factor = 1
            to_side_factor = -1
        return from_side_factor, to_side_factor
    def update_location(self, loc, increment):
        loc = [0 for i in range(self.ndim)]


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
        return ex

class LinearExtrapolateBoundaryConditionBlock(BoundaryConditionBase):
    """
    Applies 1st order linear extrapolation boundary condition.
    """
    def __init__(self, boundary_direction, side, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'LinearExtrapolate'
        return
    def halos(self):
        return True
    def apply(self, arrays, boundary_direction, side, block):
        halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
        from_side_factor, to_side_factor = self.set_side_factor(side)
        indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[boundary_direction][side]) + 1)]
        grid_vars = [GridVariable('%sG' % str(i)) for i in flatten(arrays)]
        boundary_values = [Eq(x,y) for x,y in zip(grid_vars, flatten(arrays))]
        kernel.add_equation(boundary_values)
        arrs = flatten(arrays)
        rhs_values = [2.0*i - j for i,j in zip(grid_vars, arrs)]
        equations = self.create_boundary_equations(boundary_direction, arrs, rhs_values, indices)
        kernel.add_equation(equations)
        kernel.update_block_datasets(block)
        return kernel

class DirichletBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
    """Applies constant value Dirichlet boundary condition."""
    def __init__(self, boundary_direction, side , equations, scheme= None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'Dirichlet'
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return
    def halos(self):
        return True
    def apply(self, arrays, boundary_direction, side, block):
        halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
        # Add halos across the boundary side only
        kernel.halo_ranges[boundary_direction][side] = block.boundary_halos[boundary_direction][side]
        kernel.add_equation(self.equations)
        kernel.update_block_datasets(block)
        return kernel

class SymmetryBoundaryConditionBlock(BoundaryConditionBase):
    def __init__(self, boundary_direction, side, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'Symmetry'
        return

    def apply(self, arrays, boundary_direction, side, block):
        # metric = zeros(block.ndim)
        # for i in range(block.ndim):
        #     metric[i,i] = DataObject('D%d%d' % (i,i))
        # for row in range(block.ndim):
        #     total = []
        #     expr = metric[row, :]
        #     for item in expr:
        #         total += [item**2]
        #     metric[row, :] = metric[row, :] / sqrt(sum(total))
        metric = eye(block.ndim)
        lhs_eqns = flatten(arrays)
        rhs_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                transformed_vector = metric*Matrix(ar)
                transformed_vector[boundary_direction] = -1*transformed_vector[boundary_direction]
                rhs_eqns += flatten(transformed_vector)
            else:
                rhs_eqns += [ar]

        halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
        from_side_factor, to_side_factor = self.set_side_factor(side)

        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[boundary_direction][side]) + 1)]
        final_equations = self.create_boundary_equations(boundary_direction, lhs_eqns, rhs_eqns, transfer_indices)
        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)
        return kernel

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
            f_matrix = f_matrix.transpose()[:,0:4]
        else:
            raise NotImplementedError("Side must be 0 or 1")
        return f_matrix

    def weight_function_points(self, func_points, direction, order, block, side, char_BC=False):
        if order == 1:
            h = (block.deltas[direction])**(-1)
            if side == 1:
                h = -S.One*h # Modify the first derivatives for side ==1
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
        derivatives = self.weight_function_points(fn_pts, direction, order, block, side)
        idx = block.grid_indexes[direction]
        if side == 0:
            mul_factor = 1
            start = block.ranges[direction][side]
        else:
            mul_factor = -1
            start = block.ranges[direction][side] - 1
        ecs = []
        ranges = []
        for no,d in enumerate(derivatives):
            loc = start + mul_factor*no
            ranges += [loc]
            ecs += [ExprCondPair(d, Equality(idx,loc))]
        if side != 0:
            ranges = list(reversed(ranges))
        return ecs, ranges

    def expr_cond_pair_kernel(self, fn, direction, side, order,block):
        ker = Kernel(block)
        ker.add_equation(expr)
        ker.set_computation_name("Carpenter scheme %s "%(fn))
        ker.set_grid_range(block)
        # modify the range to the number of points 
        ecs, ranges = self,expr_cond_pairs(fn, direction, side, order, block)
        ker.ranges[direction] = [ranges[0], ranges[-1]]
        return ker, ecs

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
        # for i in range(bc4.shape[0]):
        #     for j in range(bc4.shape[1]):
        #         bc4[i,j] = nsimplify(bc4[i,j])
        return bc4

class AdiabaticWallBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
    def __init__(self, boundary_direction, side, equations=None, scheme = None,plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'AdiabaticWall'
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return
    def apply(self, arrays, boundary_direction, side, block):
        halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
        wall_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                rhs = [0 for i in range(len(ar))]
                wall_eqns += [Eq(x,y) for (x,y) in zip(ar, rhs)]
        kernel.add_equation(wall_eqns)
        from_side_factor, to_side_factor = self.set_side_factor(side)
        rhs_eqns = []
        lhs_eqns = flatten(arrays)
        for ar in arrays:
            if isinstance(ar, list):
                transformed_vector = -1*Matrix(ar)
                rhs_eqns += flatten(transformed_vector)
            else:
                rhs_eqns += [ar]
        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[boundary_direction][side]) + 1)]
        final_equations = self.create_boundary_equations(boundary_direction, lhs_eqns, rhs_eqns, transfer_indices)
        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)
        return kernel

def apply_modify_derivative(order, fn, bcs, block, value):
    return

class IsothermalWallBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
    def __init__(self, boundary_direction, side, equations, const_dict, scheme = None,plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'IsothermalWall'
        self.const_dict = const_dict
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return
    def apply(self, arrays, boundary_direction, side, block):

        kernel = Kernel(block, computation_name="Isothermal wall boundary dir%d side%d" % (boundary_direction, side))
        kernel.set_boundary_plane_range(block, boundary_direction, side)
        halos = kernel.get_plane_halos(block)
        # Add the halos to the kernel in directions not equal to boundary direction
        for i in [x for x in range(block.ndim) if x != boundary_direction]:
            kernel.halo_ranges[i][0] = block.boundary_halos[i][0]
            kernel.halo_ranges[i][1] = block.boundary_halos[i][1]
        base = block.ranges[boundary_direction][side]

        metric = eye(block.ndim)
        # Set equations for the wall condition and halo points
        lhs_eqns = flatten(arrays)

        # Set wall conditions:
        wall_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                rhs = [0 for i in range(len(ar))]
                wall_eqns += [Eq(x,y) for (x,y) in zip(ar, rhs)]
        rhoE_wall = self.equations[:]
        wall_eqns += rhoE_wall
        kernel.add_equation(wall_eqns)

        n_halos = abs(halos[boundary_direction][side])
        p0 = GridVariable('p0')
        gama, Minf = ConstantObject('gama'), ConstantObject('Minf')
        vel_comps = [i**2 for i in flatten(arrays[1])]
        p0_rhs = (gama-1.0)*(lhs_eqns[-1] - 0.5*sum(vel_comps)/lhs_eqns[0])
        kernel.add_equation(Eq(p0, p0_rhs))

        loc = list(lhs_eqns[0].indices)

        for i in range(1, n_halos+1):
            new_loc = loc[:]
            new_loc[boundary_direction] += i
            T = (lhs_eqns[-1] - 0.5*(sum(vel_comps)/lhs_eqns[0]))*(gama*(gama-1.0)*Minf**2) / lhs_eqns[0]
            T = self.convert_dataset_base_expr_to_datasets(T,new_loc)
            kernel.add_equation(Eq(GridVariable('T%d' % i), T))

        if side == 0:
            from_side_factor = -1
            to_side_factor = 1
        elif side == 1:
            from_side_factor = 1
            to_side_factor = -1

        new_loc[boundary_direction] += to_side_factor
        rhs_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                transformed_vector = metric*Matrix(ar)
                transformed_vector = -1*transformed_vector
                rhs_eqns += flatten(transformed_vector)
            else:
                rhs_eqns += [ar]

        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, n_halos + 1)]

        final_equations = []
        for i, index in enumerate(transfer_indices):
            array_equations = []
            loc_lhs, loc_rhs = loc[:], loc[:]
            loc_lhs[boundary_direction] += index[0]
            loc_rhs[boundary_direction] += index[1]
            # Set rho RHS
            rhs_eqns[0] = p0*gama*Minf**2 / GridVariable('T%d' % (i+1))
            # Set rhoE RHS
            rhs_eqns[-1] = p0/(gama-1.0) + 0.5*sum(vel_comps)/lhs_eqns[0]
            for left, right in zip(lhs_eqns, rhs_eqns):
                left = self.convert_dataset_base_expr_to_datasets(left, loc_lhs)
                right = self.convert_dataset_base_expr_to_datasets(right, loc_rhs)
                array_equations += [Eq(left, right, evaluate=False)]
            final_equations += array_equations

        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)
        return kernel

class InletTransferBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
    """This is boundary condition should not be used until the user knows what he is doing. This is used for testing OpenSBLI
    """
    def __init__(self, boundary_direction, side, equations = None, scheme = None,plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'InletTransfer'
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        if side != 0:
            raise ValueError("Only implemented this BC for inlet side 0.")
        return

    def apply(self, arrays, boundary_direction, side, block):
        halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
        cons_vars = flatten(arrays)
        equations = self.create_boundary_equations(boundary_direction, cons_vars, cons_vars, [(0,-1)])
        kernel.add_equation(equations)
        kernel.update_block_datasets(block)
        return kernel


class InletPressureExtrapolateBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
    def __init__(self, boundary_direction, side, equations = None, scheme = None,plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'InletPressureExtrapolate'
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        if side != 0:
            raise ValueError("Only implemented this BC for inlet side 0.")
        return

    def apply(self, arrays, boundary_direction, side, block):
        halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
        crs = block.list_of_equation_classes[1].equations
        cr_dict = {}
        gama = ConstantObject('gama')
        for eqn in crs:
            cr_dict[str(eqn.lhs)] = eqn.rhs
        loc = arrays[0].indices
        # Evaluation of pressure, density, speed of sound on the boundary
        pb, rhob, ab = GridVariable('pb'), GridVariable('rhob'), GridVariable('ab')
        grid_vels = [GridVariable('ub%d' % i) for i in range(len(arrays[1]))]
        grid_vels_sq = [i**2 for i in grid_vels]
        eqns = [Eq(rhob, arrays[0])]
        eqns += [Eq(grid_vels[i], Abs(cr_dict['u%d_B%d' % (i, block.blocknumber)])) for i in range(len(arrays[1]))]
        eqns += [Eq(pb,(gama-1)*(flatten(arrays)[-1] - 0.5*rhob*sum(flatten(grid_vels_sq))))]
        eqns += [Eq(ab,(gama*pb/rhob)**0.5)]
        kernel.add_equation(eqns)
        arrs = flatten(arrays)
        locations = [-1, 0]
        vel = grid_vels[boundary_direction]
        # Conditions to be set at the boundary
        for lhs in flatten(arrays):
            ecs = []
            # ecs = [ExprCondPair(0, True)]
            rhs_values = [increment_dataset(lhs, boundary_direction, value) for value in locations]
            ecs += [ExprCondPair(rhs_values[0], GreaterThan(vel, ab))]
            ecs += [ExprCondPair(rhs_values[1], True)]
            kernel.add_equation(Eq(lhs, Piecewise(*ecs, **{'evaluate':False})))
        # Conditions set in the halos
        rhoE = arrays[-1]
        rhs_rhoE = pb/(gama-1.0) + 0.5*rhob*sum(grid_vels_sq)
        locations = [-i-1 for i in range(abs(halos[0][0]))]
        lhs_arrays = [increment_dataset(rhoE, boundary_direction, value) for value in locations]

        for i, lhs in enumerate(lhs_arrays):
            ecs = []
            # ecs = [ExprCondPair(0, True)]
            ecs += [ExprCondPair(lhs, GreaterThan(vel, ab))] # lhs == rhs
            ecs += [ExprCondPair(rhs_rhoE, True)]
            kernel.add_equation(Eq(lhs, Piecewise(*ecs, **{'evaluate':False})))
        pprint(kernel.equations)
        kernel.update_block_datasets(block)
        return kernel


# class OutletPressureExtrapolateBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
#     """This is boundary condition should not be used until the user knows what he is doing. This is used for testing OpenSBLI
#     """
#     def __init__(self, boundary_direction, side, equations = None, scheme = None,plane=True):
#         BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
#         self.bc_name = 'OutletPressureExtrapolate'
#         self.equations = equations
#         if not scheme:
#             self.modification_scheme = Carpenter()
#         else:
#             self.modification_scheme = scheme
#         if side != 1:
#             raise ValueError("Only implemented this BC for inlet side 1.")
#         return

#     def apply(self, arrays, boundary_direction, side, block):
#         halos, kernel = self.generate_boundary_kernel(boundary_direction, side, block, self.bc_name)
#         cons_vars = flatten(arrays)
#         equations = self.create_boundary_equations(boundary_direction, cons_vars, cons_vars, [(0,-1)])
#         kernel.add_equation(equations)
#         kernel.update_block_datasets(block)
#         return kernel

# class CharacteristicBoundaryConditionBlock(ModifyCentralDerivative, BoundaryConditionBase):
#     def __init__(self, boundary_direction, side, scheme = None,plane=True, Eigensystem = None):
#         BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
#         if not Eigensystem:
#             raise ValueError("Needs Eigen system")
#         self.Eigensystem = Eigensystem
#         return

#     def apply(self, arrays, boundary_direction, side, block):
#         # Remove the store in weno opensbli
#         # Things in """ """ should be removed
#         """
#         self.ndim = block.ndim
#         print "in charBC"
#         pprint(self.ndim)


#         self.ev_dict = self.CS.ev_store
#         self.LEV_dict = self.CS.LEV_store
#         self.REV_dict = self.CS.REV_store

#         pprint(self.LEV_dict)
#         pprint(self.REV_dict)
#         exit()"""
#         """Explanation
#         Steps for performing characteristic boundary conditions are

#         1. Evaluate LEV as GridVariable
#         2. Evaluate lambda as grid variables
#         I think 1 and 2 you can do it
#         3. REV should be evaluated as piecewise

#         """
#         self.create_REV_conditions()
#         return

#     def create_REV_conditions(self):

#         # create REV symbolic matrix, similar to WENO
#         suffix_name = 'ChBC'
#         # Before this we need to evaluate char_bc_u, and soon
#         # convert the REV into grid variables
#         rev_in_grid =self.Eigensystem.convert_matrix_to_grid_variable(self.Eigensystem.right_eigen_vector[self.direction] ,  suffix_name)
#         # These are used as the LHS symbols for the equations
#         # rev = self.Eigensystem.generate_grid_variable_REV(self.direction, suffix_name)
#         # ev = (self.Eigensystem.eigen_value[self.direction])
#         # rev = (self.Eigensystem.right_eigen_vector[self.direction])
#         # lev = (self.Eigensystem.left_eigen_vector[self.direction])
#         # self.dQ = Matrix(symbols('dQ0:%d' % 4))
#         # v = Matrix([symbols('e0:4')])
#         # pprint(v)
#         # final = rev*ev*lev*self.dQ
#         # final4 = (final[3,:].subs(EinsteinTerm('rho'), EinsteinTerm('den')))
#         # pprint(final4[0])
#         # from opensbli.core.codegeneration import ccode
#         # print ccode(final4[0],  settings = {'rational':True})
#         # exit()
#         # pprint(rev)
#         # pprint(rev[3,:])
#         # Create the RHS, making a piecewise function, you can create an empty piecewise function but that fails while printing
#         rhs_rev = zeros(*rev.shape)
#         for i in range(rhs_rev.shape[0]):
#             for j in range(rhs_rev.shape[1]):
#                 rhs_rev[i,j] = Piecewise((rev[i,j], True))
#         # Create eigen values as local variables if required
#         ev = self.Eigensystem.generate_grid_variable_ev(self.direction, suffix_name)
#         # Create expression condition pairs, for example  let's say EV[2] <0 then zero out characteristic
#         # zero out row, CHECK this for matrix indices WARNING
#         # Making up some pairs
#         #p = (Piecewise())
#         for s in range(rhs_rev.shape[0]):
#             condition_par = list(rhs_rev[s,2].atoms(ExprCondPair)) + [ExprCondPair(0, ev[s,s] <0)]
#             rhs_rev[s,2] = Piecewise(*condition_par)
#             condition_par = list(rhs_rev[s,0].atoms(ExprCondPair)) + [ExprCondPair(1, ev[s,s] >0)]
#             rhs_rev[s,1] = Piecewise(*condition_par)
#         pprint(rhs_rev)


#         return



#     def create_conditionals(self, conditions):
#         # Get velocity to check in conditional
#         ## CHECK IF NEED ABS velocity
#         a = Symbol('a')
#         u_norm = self.ev_sym[0,0]
#         dQ = []
#         for i, REV in enumerate(conditions):
#             equations = [Eq(Symbol('dW%d' % j), eq) for j, eq in enumerate(REV*self.dW_store)]
#             dQ.append(equations)
#             pprint(dQ[-1])

#         # Side 0 inlet
#         conds_0 = Piecewise((dQ[0], And(u_norm < 0, a < u_norm)),(dQ[1], And(u_norm >= 0, a < u_norm)), (dQ[2], And(u_norm >= 0, a > u_norm)), (0 , True))
#         # Side 1 outlet
#         conds_1 = Piecewise((dQ[3], And(u_norm >= 0, a < u_norm)),(dQ[4], And(u_norm < 0, a > u_norm)), (dQ[5], And(u_norm < 0, a < u_norm)), (0 , True))

#         # Side 0 inlet
#         conds_0 = Piecewise(('condition_%d' % 0, And(u_norm < 0, a < u_norm)),('condition_%d' % 1, And(u_norm >= 0, a < u_norm)), ('condition_%d' % 2, And(u_norm >= 0, a > u_norm)), (0 , True))
#         # Side 1 outlet
#         conds_1 = Piecewise(('condition_%d' % 3, And(u_norm >= 0, a < u_norm)),('condition_%d' % 4, And(u_norm < 0, a > u_norm)), ('condition_%d' % 5, And(u_norm < 0, a < u_norm)), (0 , True))
#         pprint(conds_0)
#         pprint(conds_1)
#         return

#     def generate_derivatives(self, side, order):
#         t = self.euler_eq.time_derivative
#         Q = self.euler_eq.vector_notation.get(t)
#         dQ = zeros(self.n_ev,1)
#         wrt = self.direction

#         for i, eq in enumerate(Q):
#             func_points = self.Carpenter.generate_function_points(eq, wrt, side, char_BC=True)
#             dQ[i] = self.Carpenter.compute_derivatives(func_points, wrt, order, char_BC=True)
#         return


#     def generate_symbolic_arrays(self, direction):
#         # Create symbolic evs
#         # self.symbolic_eigen_values("ev")
#         # Store ev values and placeholders
#         self.ev = self.eigen_value[direction]
#         self.ev_sym = self.eigenvalues_symbolic
#         pprint(self.ev)
#         pprint(self.ev_sym)
#         exit()
#         self.ev_values = self.eigen_value_evaluation_eq(self.ev, "ev")
#         self.dQ = Matrix(symbols('dQ0:%d' % self.n_ev))
#         # Store LEV values and placeholders
#         self.LEV = self.left_eigen_vector[direction]
#         self.LEV_sym = self.left_eigen_vector_symbolic
#         self.LEV_values = self.left_eigen_vectors_evaluation_eq(self.LEV)
#         # Store REV values and placeholders
#         self.REV = self.right_eigen_vector[direction]
#         self.REV_sym = self.right_eigen_vector_symbolic
#         self.REV_values = self.right_eigen_vectors_evaluation_eq(self.REV)

#         dW = zeros(self.n_ev,1)
#         for i in range(self.n_ev):
#             dW[i] = self.ev_sym[i,i]*self.LEV_sym[i,:]*self.dQ
#         self.dW_store = dW
#         return

#     def zero_out_characteristics(self):
#         n_ev = self.n_ev
#         # Left cases:
#         # u_neg_lt_a: u negative, subsonic, u+a zeroed out.
#         REV = self.REV_sym.copy()
#         REV.row_del(n_ev-2)
#         m_u_neg_lt_a = REV.row_insert(n_ev-2, zeros(1, n_ev))
#         # u_pos_lt_a: u positive, subsonic, u, u, u, u+a zeroed out.
#         REV = self.REV_sym.copy()
#         m_u_pos_lt_a = zeros(n_ev-1, n_ev)
#         m_u_pos_lt_a = m_u_pos_lt_a.row_insert(n_ev-1, REV[n_ev-1, :])
#         # u_pos_gt_a: u positive, supersonic, u, u, u, u+a, u-a zeroed out.
#         m_u_pos_gt_a = zeros(n_ev, n_ev)

#         # Right cases:
#         # u_pos_lt_a: u positive, subsonic, u-a zeroed out.
#         REV = self.REV_sym.copy()
#         REV.row_del(n_ev-1)
#         p_u_pos_lt_a = REV.row_insert(n_ev-1, zeros(1,n_ev))
#         # u_neg_gt_a: u negative, supersonic, u, u, u, u+a, u-a zeroed out.
#         p_u_neg_gt_a = zeros(n_ev, n_ev)
#         # u_neg_lt_a : u negative, subsonic, u, u, u, u+a zeroed out.
#         REV = self.REV_sym.copy()
#         p_u_neg_lt_a = zeros(n_ev-1, n_ev)
#         p_u_neg_lt_a = p_u_neg_lt_a.row_insert(n_ev-2, REV[n_ev-1, :])
#         conditions = [m_u_neg_lt_a, m_u_pos_lt_a, m_u_pos_gt_a, p_u_pos_lt_a, p_u_neg_gt_a, p_u_neg_lt_a]

#         return conditions


#     def return_to_conservative(self):
#         self.Q = self.REV*self.W
#         pprint(self.Q)
#         pprint(self.REV)
#         return

#     def test_eigvalues(self, ev, ndim):
#         if ndim == 1:
#             subs_list = {Symbol('a'): 3.5 ,Symbol('u0'): -0.5}
#         elif ndim == 2:
#             subs_list = {Symbol('a'): 3.5, Symbol('u0'): -0.5, Symbol('u1'): -0.25}
#         elif ndim == 3:
#             subs_list = {Symbol('a'): 3.5, Symbol('u0'): -0.5, Symbol('u1'): -0.25, Symbol('u2'): -0.05}

#         g = lambda x:x.subs(subs_list, evaluate=False)
#         return ev.applyfunc(g)
#     def evaluate_derivatives(self, direction):
#         t = self.euler_eq.time_derivative
#         Q = Matrix(self.euler_eq.vector_notation.get(t))
#         F = Matrix(self.euler_eq.vector_notation.get(direction))
#         dQ = Q.diff(direction)
#         dF = F.diff(direction)
#         pprint(dQ)
#         pprint(dF)
#         return