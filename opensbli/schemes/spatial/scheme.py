from sympy.calculus import finite_diff_weights
from sympy import postorder_traversal, Function, flatten, Rational, Idx, S
from sympy.core import Add, Mul
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable, DataSet
from opensbli.core.opensblifunctions import CentralDerivative
from opensbli.equation_types.opensbliequations import OpenSBLIEq, SimulationEquations
from opensbli.core.kernel import Kernel
from opensbli.utilities.helperfunctions import increasing_order
from opensbli.core.datatypes import Int


class Scheme(object):

    """ A numerical discretisation scheme. """

    def __init__(self, name, order):
        """ Initialise the scheme.

        :arg str name: The name of the scheme.
        :arg int order: The order of the scheme.
        """

        self.name = name
        self.order = order
        return


class CentralHalos(object):
    def __init__(self, order):
        # Check for the boundary types in the blocks and set the halo points
        # self.halos = [[-scheme.order, scheme.order] for dim in range(block.ndim)]
        # self.halos = [-5, 5]
        self.halos = [-order/2, order/2]
        return

    def get_halos(self, side):
        return self.halos[side]

    def __str__(self):
        return "CentralHalos"


class CentralHalos_defdec(object):
    def __init__(self):
        # Check for the boundary types in the blocks and set the halo points
        # self.halos = [[-scheme.order, scheme.order] for dim in range(block.ndim)]
        self.halos = [-5, 5]
        return

    def get_halos(self, side):
        return self.halos[side]

    def __str__(self):
        return "CentralHalos"


class Central(Scheme):

    """ Spatial discretisation scheme using central differences. """

    def __init__(self, order):
        """ Set up the scheme.

        :arg int order: The order of accuracy of the scheme.
        """
        Scheme.__init__(self, "CentralDerivative", order)
        print "A Central scheme of order %d is being used" % order
        self.schemetype = "Spatial"
        # Points for the spatial scheme
        self.points = list(i for i in range(-order/2, order/2+1))
        # Set max_order to 2 currently
        max_order = 2
        # self._generate_derivative_weights(max_order)
        self.required_constituent_relations = {}
        self.halotype = CentralHalos(order)
        return

    def _generate_weights(self, direction, order, block):
        """ Descritises only the homogeneous derivatives of any order or
        first derivatives"""
        self.diffpoints = [i for i in self.points]
        weights = finite_diff_weights(order, self.diffpoints, 0)
        return weights[order][-1]

    def add_required_database(self, dbases):
        self.required_database += flatten(list(dbases))
        return

    @property
    def scheme_required_databases(self):
        return set(self.required_database)

    def update_works(self, to_descritse, block):

        return

    def set_halos(self, block):
        for direction in range(block.ndim):
            block.set_block_boundary_halos(direction, 0, self.halotype)
            block.set_block_boundary_halos(direction, 1, self.halotype)
        return

    def discretise(self, type_of_eq, block):
        """
        This is main calling function from opensbli equations.spatial_discretisation, which is called from block.discretise
        Do here the following
        a. Find all the Function of type(self.name)
        b. Add all the DataSetBases required to the type_of_eq
        c. Create Equations for the evaluation ro create Kernels of each function depending on the grid/block
        control parameters
        d. Set the range of evaluation of the DataSetBases
        e. Update the Descritised equations in type_of_eq by substituting the equations with respective
        work array or discretised formula
        """
        # Check if it is similar to compressible Navier stokes equations
        # if type_of_eq.
        self.set_halos(block)
        if isinstance(type_of_eq, SimulationEquations):
            """ Simulation equations are always solved as sbli_rhs_discretisation as of now"""
            # if block.sbli_rhs_discretisation:
            self.sbli_rhs_discretisation(type_of_eq, block)
            return self.required_constituent_relations
        else:
            block.store_work_index # Store work
            local_kernels, discretised_eq = self.general_discretisation(type_of_eq.equations, block, name=type_of_eq.__class__.__name__)
            block.reset_work_to_stored # Reset
            if discretised_eq:
                for ker in local_kernels:
                    eval_ker = local_kernels[ker]
                    type_of_eq.Kernels += [eval_ker]

                discretisation_kernel = Kernel(block, computation_name="%s evaluation" % type_of_eq.__class__.__name__)
                discretisation_kernel.set_grid_range(block)
                for eq in discretised_eq:
                    discretisation_kernel.add_equation(eq)
                discretisation_kernel.update_block_datasets(block)
                type_of_eq.Kernels += [discretisation_kernel]
                return self.required_constituent_relations
            else:
                pass
            return self.required_constituent_relations

    def get_local_function(self, list_of_components):
        CD_fns = []
        for c in list_of_components:
            CD_fns += list(c.atoms(CentralDerivative))
        CD_fns = list(set(CD_fns))
        return CD_fns

    def group_by_direction(self, eqs):
        all_CDS = []
        for eq in eqs:
            all_CDS += list(eq.atoms(CentralDerivative))
        all_CDS = list(set(all_CDS))
        grouped = {}
        for cd in all_CDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped

    def update_range_of_constituent_relations(self, CD, block):
        direction = CD.get_direction[0]

        if CD.required_datasets:
            for v in CD.required_datasets:
                if v in self.required_constituent_relations.keys():
                    self.required_constituent_relations[v].set_halo_range(direction, 0, self.halotype)
                    self.required_constituent_relations[v].set_halo_range(direction, 1, self.halotype)
                else:
                    self.required_constituent_relations[v] = Kernel(block, computation_name="CR%s" % v)
                    self.required_constituent_relations[v].set_grid_range(block)
                    self.required_constituent_relations[v].set_halo_range(direction, 0, self.halotype)
                    self.required_constituent_relations[v].set_halo_range(direction, 1, self.halotype)
        return

    def check_constituent_relations(self, block, list_of_eq):
        """ Checks all the datasets in equations provided are evaluated in constituent relations."""
        arrays = []
        for eq in flatten(list_of_eq):
            arrays += list(eq.atoms(DataSet))
        arrays = set(arrays)
        undefined = arrays.difference(self.required_constituent_relations.keys())

        for dset in undefined:
            self.required_constituent_relations[dset] = Kernel(block, computation_name="CR%s" % dset)
            self.required_constituent_relations[dset].set_grid_range(block)
        return

    def sbli_rhs_discretisation(self, type_of_eq, block):
        """
        This is the discretisation similar to the Southampton SBLI RHS, without
        entropy splitting. It is upto the user to give the equations in Conservative
        formulation
        """
        equations = flatten(type_of_eq.equations)
        residual_arrays = [eq.residual for eq in equations]
        equations = [e._sanitise_equation for e in equations]
        classify_parameter = ConstantObject("Re")
        self.required_constituent_relations = {}
        viscous, convective = self.classify_equations_on_parameter(equations, classify_parameter)
        kernels = []
        convective_grouped = self.group_by_direction(convective)
        if convective_grouped:
            # Create equations for evaluation of derivatives
            for key, value in convective_grouped.iteritems():
                for v in value:
                    v.update_work(block)
            local_evaluations_group = {}
            function_expressions_group = {}
            # Process the convective derivatives, this requires grouping of equations
            subs_conv = {}
            for key, value in convective_grouped.iteritems():
                local_evaluations_group[key] = []
                ev_ker = Kernel(block)
                ev_ker.set_computation_name("Convective terms group %d" % key)
                ev_ker.set_grid_range(block)
                block.store_work_index
                local = []
                for v in value:
                    self.update_range_of_constituent_relations(v, block)
                    if len(v.required_datasets) > 1:
                        wk = block.work_array()
                        block.increase_work_index
                        expr = OpenSBLIEq(wk, v.args[0])
                        ev_ker.add_equation(expr)
                        ev_ker.set_halo_range(key, 0, self.halotype)
                        ev_ker.set_halo_range(key, 1, self.halotype)
                        v1 = v.subs(v.args[0], wk)
                    else:
                        v1 = v
                    expr = OpenSBLIEq(v.work, v1._discretise_derivative(self, block))
                    ker = Kernel(block)
                    ker.add_equation(expr)
                    ker.set_computation_name("Convective %s " % (v))
                    ker.set_grid_range(block)
                    local += [ker]
                    subs_conv[v] = v.work
                if ev_ker.equations:
                    local_evaluations_group[key] += [ev_ker]
                function_expressions_group[key] = local
                block.reset_work_to_stored
            # Convective evaluation
            for key, value in local_evaluations_group.iteritems():
                kernels += value + function_expressions_group[key]
            # Create convective residual

            convective_descritised = convective[:]
            self.check_constituent_relations(block, convective_descritised)

            for no, c in enumerate(convective_descritised):
                convective_descritised[no] = convective_descritised[no].subs(subs_conv)
            conv_residual_kernel = self.create_residual_kernel(residual_arrays, convective_descritised, block)
            conv_residual_kernel.set_computation_name("Convective residual ")
            kernels += [conv_residual_kernel]
        # reset the work index of blocks
        block.reset_work_index
        # Discretise the viscous fluxes. This is straight forward as we need not modify anything
        viscous_kernels, viscous_discretised = self.general_discretisation(viscous, block, name="Viscous")
        self.check_constituent_relations(block, viscous)
        if viscous_kernels:
            for ker in sorted(viscous_kernels, cmp=increasing_order):
                eval_ker = viscous_kernels[ker]
                kernels += [eval_ker]
        if viscous_discretised:
            visc_residual_kernel = self.create_residual_kernel(residual_arrays, viscous_discretised, block)
            visc_residual_kernel.set_computation_name("Viscous residual")
            kernels += [visc_residual_kernel]
        # Add the kernels to the solutions
        type_of_eq.Kernels += kernels
        block.reset_work_index
        return

    def set_halo_range_kernel(self, kernel, direction, sides=None):
        if not sides:
            kernel.set_halo_range(direction, 0, self.halotype)
            kernel.set_halo_range(direction, 1, self.halotype)
            return kernel
        else:
            raise NotImplementedError("")

    def create_residual_kernel(self, residual_arrays, discretised_eq, block):
        if len(residual_arrays) != len(discretised_eq):
            raise ValueError("")
        residue_kernel = Kernel(block)
        for no, array in enumerate(residual_arrays):
            expr = OpenSBLIEq(array, array+discretised_eq[no])
            residue_kernel.add_equation(expr)
        residue_kernel.set_grid_range(block)
        return residue_kernel

    def general_discretisation(self, equations, block, name=None):
        """
        This discretises the central derivatives, without any special treatment to
        group the derivatives or any thing
        """
        descritised_equations = flatten(equations)[:]
        cds = self.get_local_function(flatten(equations))
        if cds:
            local_kernels = {}
            if block.store_derivatives:
                for der in cds:
                    der.update_work(block)
                    ker = Kernel(block)
                    if name:
                        ker.set_computation_name("%s %s " % (name, der))
                    local_kernels[der] = ker  # Reverted back
            # create a dictionary of works and kernels
            work_arry_subs = {}
            for der in cds:
                self.update_range_of_constituent_relations(der, block)
                expr, local_kernels = self.traverse(der, local_kernels, block)
                expr_discretised = OpenSBLIEq(der.work, expr._discretise_derivative(self, block))
                work_arry_subs[expr] = der.work
                local_kernels[der].add_equation(expr_discretised)
                local_kernels[der].set_grid_range(block)
            for no, c in enumerate(descritised_equations):
                descritised_equations[no] = descritised_equations[no].subs(work_arry_subs)
            return local_kernels, descritised_equations
        else:
            return None, None

    def classify_equations_on_parameter(self, equations, parameter):
        containing_terms = [S.Zero for eq in equations]
        other_terms = [S.Zero for eq in equations]
        for number, eq in enumerate(equations):
            if isinstance(eq.rhs, Add):
                for expr in eq.rhs.args:
                    if expr.has(parameter):
                        containing_terms[number] = Add(containing_terms[number], expr)
                    else:
                        other_terms[number] = Add(other_terms[number], expr)
            elif isinstance(eq.rhs, Mul):
                expr = eq.rhs
                if expr.has(parameter):
                    containing_terms[number] = Add(containing_terms[number], expr)
                else:
                    other_terms[number] = Add(other_terms[number], expr)
        # Zero out any other derivatives in the containing terms and other terms
        for no, eq in enumerate(other_terms):
            fns = [fn for fn in eq.atoms(Function) if not isinstance(fn, CentralDerivative)]
            substitutions = dict(zip(fns, [0]*len(fns)))
            other_terms[no] = other_terms[no].subs(substitutions)

        for no, eq in enumerate(containing_terms):
            fns = [fn for fn in eq.atoms(Function) if not isinstance(fn, CentralDerivative)]
            substitutions = dict(zip(fns, [0]*len(fns)))
            containing_terms[no] = containing_terms[no].subs(substitutions)
        return containing_terms, other_terms

    def traverse(self, CD, kernel_dictionary, block):
        expr = CD.copy()
        inner_cds = []
        # if CD.args[0].atoms(CentralDerivative):
        pot = postorder_traversal(CD)
        inner_cds = []
        for p in pot:
            if isinstance(p, CentralDerivative):
                inner_cds += [p]
            else:
                continue
        # Contains inner derivatives
        if len(inner_cds) > 1:
            for np, cd in enumerate(inner_cds[:-1]):
                if cd.is_store and cd.work:
                    cd.is_used(True)
                    expr = expr.subs(cd, cd.work)
                    # update the kernel ranges for inner cds
                    dires = inner_cds[np+1].get_direction
                    for direction in dires:
                        kernel_dictionary[cd].set_halo_range(direction, 0, self.halotype)
                        kernel_dictionary[cd].set_halo_range(direction, 1, self.halotype)
                elif cd.is_store:
                    raise ValueError("NOT IMPLEMENTED THIS")
                elif not cd.is_store:  # Donot store derivatives
                    raise ValueError("This dependency sgould be validated for Carpenter BC")
                    expr = expr.subs(cd, cd._discretise_derivative(self, block))
                else:
                    raise ValueError("Could not classify this")
        return expr, kernel_dictionary


class TemporalSolution(object):
    def __init__(self):
        self.start_kernels = []
        self.end_kernels = []
        self.kernels = []
        return
