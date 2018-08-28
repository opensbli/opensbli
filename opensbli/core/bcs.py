"""@brief
   @authors David J Lusher, Satya Pramod Jammy
   @contributors
   @details
"""
from sympy import flatten, zeros, Matrix, S, sqrt, Equality, nsimplify, Float, Rational, Idx, pprint
from opensbli.core.kernel import Kernel, ConstantsToDeclare
from opensbli.core.opensbliobjects import DataSet, ConstantIndexed, ConstantObject
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.datatypes import Int
from opensbli.core.grid import GridVariable
from opensbli.physical_models.ns_physics import NSphysics
from opensbli.utilities.helperfunctions import get_min_max_halo_values, increment_dataset
from sympy.functions.elementary.piecewise import ExprCondPair
from opensbli.schemes.spatial.weno import ShockCapturing
side_names = {0: 'left', 1: 'right'}
# V2: TODO, decide which boundary conditions should be one sided or not


class Exchange(object):
    pass


class ModifyCentralDerivative(object):
    """ A place holder for the boundary conditions on which the central derivative should be modified"""
    pass


class ExchangeSelf(Exchange):
    """ Defines data exchange on the same block. """

    def __init__(self, block, direction, side):
        # Range of evaluation (i.e. the grid points, including the halo points, over which the computation should be performed).
        self.computation_name = "exchange"
        self.block_number = block.blocknumber
        self.block_name = block.blockname
        self.direction = direction
        self.side = side_names[side]
        return

    @property
    def name(self):
        return "%s%d" % (self.computation_name, self.number)

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
        name = "Boundary_exchange_block_%d_direction%d_side%s" % (self.block_number, self.direction, self.side)
        return name

    def write_latex(self, latex):
        latex.write_string("This is an exchange self kernel on variables %s\\\\" % ', '.join([str(a) for a in self.transfer_arrays]))
        return

    @property
    def opsc_code(self):
        """The string for calling the boundary condition in OPSC is updated while creating
        the code for exchanging data."""
        return [self.call_name]


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
        self.boundary_types = zip(it, it)
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
        print kernel.computation_name
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


class PeriodicBC(BoundaryConditionBase):
    """ Applies an exchange periodic boundary condition.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
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
        boundary_direction, side = self.direction, self.side
        transfer_from = [d[0] for d in halos]
        transfer_to = [d[0] for d in halos]
        if side == 0:
            transfer_from[boundary_direction] = idx[boundary_direction].lower
            transfer_to[boundary_direction] = idx[boundary_direction].upper
        else:
            transfer_from[boundary_direction] = idx[boundary_direction].upper + halos[boundary_direction][0]
            transfer_to[boundary_direction] = idx[boundary_direction].lower + halos[boundary_direction][0]

        transfer_size = Matrix([i.upper + i.lower for i in idx]) + \
            Matrix([abs(dire[0]) + abs(dire[1]) for dire in halos])
        transfer_size[boundary_direction] = abs(halos[boundary_direction][side])
        return transfer_size, transfer_from, transfer_to


class DirichletBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Applies a constant value Dirichlet boundary condition.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg list equations: OpenSBLI equations to enforce on the boundary.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, equations, scheme=None, plane=True):
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

    def apply(self, arrays, block):
        direction, side = self.direction, self.side
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        # Dirichlet set on the boundary
        kernel.add_equation(self.equations)
        # Change ranges if using split BC
        if isinstance(kernel.halo_ranges, ConstantIndexed):
            # Manually set Dirichlet into the halo range of this side
            halo_object = kernel.halo_ranges
            halo_object._value[2*direction + side] = halos[direction][side]
        else:  # Not using split BC, halos should be updated
            kernel.halo_ranges[direction][side] = block.boundary_halos[direction][side]
        kernel.update_block_datasets(block)
        return kernel


class SymmetryBC(BoundaryConditionBase):
    """ Applies a symmetry condition on the boundary, normal velocity components set/evaluate to zero.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'Symmetry'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        fd_metric = block.fd_metrics
        direction, side = self.direction, self.side
        direction_metric = Matrix(fd_metric[direction, :])
        normalisation = sqrt(sum([a**2 for a in direction_metric]))
        unit_normals = direction_metric/normalisation
        lhs_eqns = flatten(arrays)
        boundary_values = []
        rhs_eqns = []
        for ar in arrays:
            if isinstance(ar, list):
                contra_variant_vector = unit_normals.dot(Matrix(ar))
                transformed_vector = Matrix(ar).T - 2.*contra_variant_vector*unit_normals
                rhs_eqns += flatten(transformed_vector)
                # Later move this to an inviscid wall boundary condition
                transformed_vector = Matrix(ar).T - 1.*contra_variant_vector*unit_normals
                boundary_values += flatten(transformed_vector)
            else:
                rhs_eqns += [ar]
                boundary_values += [ar]
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        from_side_factor, to_side_factor = self.set_side_factor()

        transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[direction][side]) + 1)]
        final_equations = self.create_boundary_equations(lhs_eqns, rhs_eqns, transfer_indices)
        transfer_indices = [tuple([0, to_side_factor])]
        final_equations += self.create_boundary_equations(lhs_eqns, boundary_values, transfer_indices)
        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)
        return kernel


class Carpenter(object):
    """ 4th order one-sided Carpenter boundary treatment (https://doi.org/10.1006/jcph.1998.6114).
    If a boundary condition is an instance of ModifyCentralDerivative,
    central derivatives are replaced at that domain boundary by the Carpenter scheme."""

    def __init__(self):
        self.bc4_coefficients = self.carp4_coefficients()
        self.bc4_2_coefficients = self.second_der_coefficients()
        return

    def function_points(self, expression, direction, side):
        """ Create the function locations for evaluation of the Carpenter
        one-sided derivatives."""
        f_matrix = zeros(6, 6)
        loc = list(list(expression.atoms(DataSet))[0].indices)
        for shift in range(6):
            func_points = []
            for index in range(6):
                new_loc = loc[:]
                new_loc[direction] += index - shift
                for dset in expression.atoms(DataSet):
                    expression = expression.replace(dset, dset.base[new_loc])
                func_points.append(expression)
            f_matrix[:, shift] = func_points
        if side == 0:
            f_matrix = f_matrix[:, 0:4]
        elif side == 1:
            f_matrix = f_matrix.transpose()[:, 0:4]
        else:
            raise NotImplementedError("Side must be 0 or 1")
        return f_matrix

    def weight_function_points(self, func_points, direction, order, block, side):
        """ Multiply the function points with their weightings to form the derviatives."""
        if order == 1:
            h = S.One  # The division of delta is now applied in central derivative to reduce the divisions
            if side == 1:
                h = -S.One*h  # Modify the first derivatives for side ==1
            weighted = zeros(4, 1)
            for i in range(4):
                weighted[i] = h*(self.bc4_coefficients[i, :]*func_points[:, i])
        elif order == 2:
            h_sq = S.One  # The division of delta**2 is now applied in central derivative to reduce the divisions
            weighted = zeros(2, 1)
            for i in range(2):
                weighted[i] = h_sq*(self.bc4_2_coefficients[i, :]*func_points[0:5, i])
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
        for no, d in enumerate(derivatives):
            loc = start + mul_factor*no
            ecs += [ExprCondPair(d, Equality(idx, loc))]
        return ecs

    def second_der_coefficients(self):
        """ Computes the finite-difference coefficients for the 2nd order one sided Carpenter wall boundary derivative.
        returns: Matrix: bc4_2: Matrix of stencil coefficients."""
        bc4_2 = Matrix([[35.0, -104.0, 114.0, -56.0, 11.0], [11.0, -20.0, 6.0, 4.0, -1.0]])/12.0
        for i in range(bc4_2.shape[0]):
            for j in range(bc4_2.shape[1]):
                bc4_2[i, j] = nsimplify(bc4_2[i, j])
        return bc4_2

    def carp4_coefficients(self):
        """ Computes the finite-difference coefficients for the 1st order one sided Carpenter wall boundary derivative.
        :returns: Matrix: bc4: Matrix of stencil coefficients."""
        R1 = -(2177.0*sqrt(295369.0)-1166427.0)/25488.0
        R2 = (66195.0*sqrt(53.0)*sqrt(5573.0)-35909375.0)/101952.0

        al4_0 = [-(216.0*R2+2160.0*R1-2125.0)/12960.0, (81.0*R2+675.0*R1+415.0)/540.0, -(72.0*R2+720.0*R1+445.0)/1440.0, -(108.0*R2+756.0*R1+421.0)/1296.0]
        al4_1 = [(81.0*R2+675.0*R1+415.0)/540.0, -(4104.0*R2+32400.0*R1+11225.0)/4320.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(216.0*R2+2160.0*R1+655.0)/4320.0]
        al4_2 = [-(72.0*R2+720.0*R1+445.0)/1440.0, (1836.0*R2+14580.0*R1+7295.0)/2160.0, -(4104.0*R2+32400.0*R1+12785.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0]
        al4_3 = [-(108.0*R2+756.0*R1+421.0)/1296.0, -(216.0*R2+2160.0*R1+655.0)/4320.0, (81.0*R2+675.0*R1+335.0)/540.0, -(216.0*R2+2160.0*R1-12085.0)/12960.0]

        al4 = Matrix([al4_0, al4_1, al4_2, al4_3])

        ar4_0 = [(-1.0)/2.0, -(864.0*R2+6480.0*R1+305.0)/4320.0, (216.0*R2+1620.0*R1+725.0)/540.0, -(864.0*R2+6480.0*R1+3335.0)/4320.0, 0.0, 0.0]
        ar4_1 = [(864.0*R2+6480.0*R1+305.0)/4320.0, 0.0, -(864.0*R2+6480.0*R1+2315.0)/1440.0, (108.0*R2+810.0*R1+415.0)/270.0, 0.0, 0.0]
        ar4_2 = [-(216.0*R2+1620.0*R1+725.0)/540.0, (864.0*R2+6480.0*R1+2315.0)/1440.0, 0.0, -(864.0*R2+6480.0*R1+785.0)/4320.0, -1.0/12.0, 0.0]
        ar4_3 = [(864.0*R2+6480.0*R1+3335.0)/4320.0, -(108.0*R2+810.0*R1+415.0)/270.0, (864.0*R2+6480.0*R1+785.0)/4320.0, 0.0, 8.0/12.0, -1.0/12.0]
        ar4 = Matrix([ar4_0, ar4_1, ar4_2, ar4_3])
        # Form inverse and convert to rational
        al4_inv = al4.inv()
        bc4 = al4_inv*ar4
        return bc4


class IsothermalWallBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Navier-Stokes specific boundary condition. Applies a no-slip viscous wall condition,
    velocity components are zero on the wall. Temperature is fixed with a prescribed wall temperature,
    given as the rhoE equation passed to this BC.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg Eq equations: Equation for conservative variable rhoE to set on the wall with constant temperature.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, equations, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'IsothermalWall'
        self.equations = equations
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        n_halos = abs(halos[self.direction][self.side])
        # Using Navier Stokes physics object, create conservative variables
        NS = NSphysics(block)
        from_side_factor, to_side_factor = self.set_side_factor()
        wall_eqns = []
        # Set wall conditions, momentum zero
        wall_eqns += [OpenSBLIEq(x, Float(S.Zero)) for x in NS.momentum()]
        # Set wall energy, density is left to be evaluated
        wall_eqns += self.equations[:]
        kernel.add_equation(wall_eqns)
        # Update halos if a shock capturing scheme is being used.
        if any(isinstance(sc, ShockCapturing) for sc in block.discretisation_schemes.values()):
            pw, gamma, Minf, Twall = GridVariable('Pwall'), NS.specific_heat_ratio(), NS.mach_number(), ConstantObject('Twall')
            # Wall pressure
            wall_eqns = [OpenSBLIEq(pw, increment_dataset(NS.pressure(relation=True, conservative=True), self.direction, 0))]
            kernel.add_equation(wall_eqns)
            # Direct copy indices over the wall
            transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, n_halos + 1)]
            # Evaluate velocities in the halos
            halo_velocities = []
            for index, momentum_comp in enumerate(NS.momentum()):
                for pair_index, index_pair in enumerate(transfer_indices):
                    momentum, density = increment_dataset(momentum_comp, self.direction, index_pair[1]), increment_dataset(NS.density(), self.direction, index_pair[1])
                    halo_velocities += [OpenSBLIEq(GridVariable('u%d%d' % (index, pair_index+1)), momentum/density)]
            kernel.add_equation(halo_velocities)
            # Evaluate the temperature one point above the wall, use it to interpolate temperature in the halos
            T_above = GridVariable('T_above')
            kernel.add_equation(OpenSBLIEq(GridVariable('T_above'), increment_dataset(NS.temperature(relation=True, conservative=True), self.direction, to_side_factor)))
            halo_densities = []
            halo_momentum = []
            halo_energy = []
            # Evaluate halo density using wall pressure and interpolated temperatures

            for i in range(1, n_halos+1):
                # Interpolate temperature into the halos and calculate halo density
                rho_halo = GridVariable('rho_halo_%d' % i)
                T_halo = GridVariable('T%d' % i)
                kernel.add_equation([OpenSBLIEq(T_halo, (i+1)*Twall - i*T_above)])
                halo_densities += [OpenSBLIEq(rho_halo, pw*gamma*Minf*Minf/T_halo)]
                # Set momentum in the halos, reverse direction
                velocities = []
                for index, momentum_comp in enumerate(NS.momentum()):
                    velocity = GridVariable('u%d%d' % (index, i))
                    halo_momentum += [OpenSBLIEq(increment_dataset(momentum_comp, self.direction, from_side_factor*i), -rho_halo*velocity)]
                    velocities += [velocity**2]
                # Set energy in halos based on wall pressure, interpolated temperature and the calculated density
                energy_rhs = pw/(gamma-1) + Rational(1, 2)*rho_halo*(sum(velocities))
                halo_energy += [OpenSBLIEq(increment_dataset(NS.total_energy(), self.direction, from_side_factor*i), energy_rhs)]
                halo_densities += [OpenSBLIEq(increment_dataset(NS.density(), self.direction, from_side_factor*i), rho_halo)]
            kernel.add_equation(halo_densities)
            kernel.add_equation(halo_momentum)
            kernel.add_equation(halo_energy)
        kernel.update_block_datasets(block)
        return kernel


class InletTransferBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Simple inlet boundary condition to copy all solution variable values from the left halos
    to the boundary plane.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'InletTransfer'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        if side != 0:
            raise ValueError("Only implemented this BC for inlet side 0.")
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        cons_vars = flatten(arrays)
        equations = self.create_boundary_equations(cons_vars, cons_vars, [(0, -1)])
        kernel.add_equation(equations)
        kernel.update_block_datasets(block)
        return kernel


class ExtrapolationBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Extrapolation boundary condition. Pass order 0 for Zeroth order extrapolation
    and order=1 for linear extrapolation.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only.
    TODO Is it modify central"""

    def __init__(self, boundary_direction, side, order, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'Extrapolation'
        # Order of the extrapolation
        self.order = order
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        cons_vars = flatten(arrays)
        n_halos = abs(halos[self.direction][self.side])
        from_side_factor, to_side_factor = self.set_side_factor()

        if self.order == 0:  # Zeroth order extrapolation
            halo_points = [from_side_factor*i for i in range(0, n_halos+1)]
            for i in halo_points:
                equations = self.create_boundary_equations(cons_vars, cons_vars, [(i, to_side_factor)])
                kernel.add_equation(equations)
        else:
            raise ValueError("Only zeroth order currently implemented.")
        kernel.update_block_datasets(block)
        return kernel


class AdiabaticWallBC(ModifyCentralDerivative, BoundaryConditionBase):
    """ Adiabatic wall condition, zero gradient dT/dn = 0 over the boundary.

    :arg int boundary_direction: Spatial direction to apply boundary condition to.
    :arg int side: Side 0 or 1 to apply the boundary condition for a given direction.
    :arg object scheme: Boundary scheme if required, defaults to Carpenter boundary treatment.
    :arg bool plane: True/False: Apply boundary condition to full range/split range only."""

    def __init__(self, boundary_direction, side, scheme=None, plane=True):
        BoundaryConditionBase.__init__(self, boundary_direction, side, plane)
        self.bc_name = 'AdiabaticWall'
        if not scheme:
            self.modification_scheme = Carpenter()
        else:
            self.modification_scheme = scheme
        return

    def apply(self, arrays, block):
        halos, kernel = self.generate_boundary_kernel(block, self.bc_name)
        wall_eqns = []
        from_side_factor, to_side_factor = self.set_side_factor()
        for ar in arrays:
            if isinstance(ar, list):  # Set velocity components to zero on the wall
                rhs = [0 for i in range(len(ar))]
                wall_eqns += [OpenSBLIEq(x, y) for (x, y) in zip(ar, rhs)]
            else:  # Take rho and rhoE from one point above the wall
                wall_eqns += [OpenSBLIEq(ar, increment_dataset(ar, self.direction, to_side_factor))]
        kernel.add_equation(wall_eqns)
        final_equations = []
        if any(isinstance(sc, ShockCapturing) for sc in block.discretisation_schemes.values()):
            rhs_eqns = []
            lhs_eqns = flatten(arrays)
            for ar in arrays:
                if isinstance(ar, list):  # Velocity components in the halos are set negative
                    transformed_vector = -1*Matrix(ar)
                    rhs_eqns += flatten(transformed_vector)
                else:  # rho and rhoE are copied symmetrically over the boundary
                    rhs_eqns += [ar]
            transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[self.direction][self.side]) + 1)]
            final_equations = self.create_boundary_equations(lhs_eqns, rhs_eqns, transfer_indices)
        kernel.add_equation(final_equations)
        kernel.update_block_datasets(block)
        return kernel
