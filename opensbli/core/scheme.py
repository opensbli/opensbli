
from sympy.calculus import finite_diff_weights
from sympy import *
from sympy.printing import *
from .opensbliobjects import ConstantObject
from .opensblifunctions import CentralDerivative
from .opensbliequations import OpenSBLIExpression
from .weno_opensbli import *
from .kernel import *
from .latex import *
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
        #self.halos = [[-scheme.order, scheme.order] for dim in range(block.ndim)]
        self.halos = [-5, 5]
        #self.halos = [-order/2, order/2]
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
        self.schemetype = "Spatial"
        # Points for the spatial scheme
        self.points = list(i for i in range(-order/2, order/2+1))
        # Set max_order to 2 currently
        max_order = 2
        #self._generate_derivative_weights(max_order)
        self.required_constituent_relations = {}
        self.halotype = CentralHalos(order)
        return

    def _generate_weights(self, direction, order, block):
        """ Descritises only the homogeneous derivatives of any order or
        first derivatives"""
        #print(self.points, direction, type(direction))
        self.diffpoints = [i*block.deltas[direction] for i in self.points]
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
            block.set_block_boundary_halos(direction,0, self.halotype )
            block.set_block_boundary_halos(direction, 1, self.halotype )
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
        #if type_of_eq.
        self.set_halos(block)
        from .opensbliequations import *
        if isinstance(type_of_eq, SimulationEquations):
            """ Simulation equations are always solved as sbli_rhs_discretisation as of now"""
            #if block.sbli_rhs_discretisation:
            self.sbli_rhs_discretisation(type_of_eq, block)
            return self.required_constituent_relations
        else:
            local_kernels, discretised_eq = self.genral_discretisation(type_of_eq.equations, block, name =  type_of_eq.__class__.__name__)
            if discretised_eq:
                for ker in local_kernels:
                    eval_ker = local_kernels[ker]
                    #eval_ker.set_computation_name("%s "%(ker))
                    #eval_ker.update_block_datasets(block)
                    type_of_eq.Kernels += [eval_ker]

                discretisation_kernel = Kernel(block, computation_name="%s evaluation"%type_of_eq.__class__.__name__)
                discretisation_kernel.set_grid_range(block)
                for eq in discretised_eq:
                    discretisation_kernel.add_equation(eq)
                discretisation_kernel.update_block_datasets(block)
                type_of_eq.Kernels +=  [discretisation_kernel]
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
                    self.required_constituent_relations[v] = Kernel(block, computation_name = "CR%s"%v)
                    self.required_constituent_relations[v].set_grid_range(block)
                    self.required_constituent_relations[v].set_halo_range(direction, 0, self.halotype)
                    self.required_constituent_relations[v].set_halo_range(direction, 1, self.halotype)
        return

    def sbli_rhs_discretisation(self, type_of_eq, block):
        """
        This is the discretisation similar to the Southampton SBLI RHS, without
        entropy splitting. It is upto the user to give the equations in Conservative
        formulation
        """
        equations =  flatten(type_of_eq.equations)
        residual_arrays = [eq.residual for eq in equations]
        equations = [e._sanitise_equation for e in equations]
        classify_parameter = ConstantObject("Re")
        self.required_constituent_relations = {}
        viscous, convective = self.classify_equations_on_parameter(equations, classify_parameter)
        convective_grouped = self.group_by_direction(convective)
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
            ev_ker.set_computation_name("Convective terms group %d"%key)
            ev_ker.set_grid_range(block)
            block.store_work_index
            local = []
            for v in value:
                self.update_range_of_constituent_relations( v, block)
                if len(v.required_datasets) >1:
                    wk = block.work_array()
                    block.increase_work_index
                    expr = Eq(wk, v.args[0])
                    ev_ker.add_equation(expr)
                    ev_ker.set_halo_range(key, 0, self.halotype)
                    ev_ker.set_halo_range(key, 1, self.halotype)
                    #expr.work = wk
                    v1 = v.subs(v.args[0], wk)
                else:
                    v1 = v
                expr = Eq(v.work, v1._discretise_derivative(self, block))
                #pprint(expr)
                ker = Kernel(block)
                ker.add_equation(expr)
                ker.set_computation_name("Convective %s "%(v))
                ker.set_grid_range(block)
                #pprint(ker.__dict__)
                #expr.work = v.work
                local += [ker]
                # create any Boundary modification kernels Reverted back
                #local += v1.apply_boundary_derivative_modification(block, self, v.work)
                subs_conv[v] = v.work
            #pprint(ev_ker.__dict__)
            if ev_ker.equations:
                local_evaluations_group[key] += [ev_ker]
            function_expressions_group[key] = local
            block.reset_work_to_stored
        # To do evaluate the convective residuals
        kernels = []
        #Convective evaluation
        for key, value in local_evaluations_group.iteritems():
            kernels += value + function_expressions_group[key]
        # Create convective residual

        convective_descritised = convective[:]

        for no, c in enumerate(convective_descritised):
            convective_descritised[no] = convective_descritised[no].subs(subs_conv)

        conv_residual_kernel = self.create_residual_kernel(residual_arrays,convective_descritised, block)
        conv_residual_kernel.set_computation_name("Convective residual ")
        kernels += [conv_residual_kernel]
        # reset the work index of blocks
        block.reset_work_index
        # Discretise the viscous fluxes. This is straight forward as we need not modify anything
        viscous_kernels, viscous_discretised = self.genral_discretisation(viscous, block, name="Viscous")
        if viscous_kernels:
            for ker in viscous_kernels:
                eval_ker = viscous_kernels[ker]
                #eval_ker.set_computation_name("Viscous %s "%(ker))
                kernels += [eval_ker]
                #pprint(['Viscous things going', ker])
                #kernels += ker.apply_boundary_derivative_modification(block, self, ker.work)
        if viscous_discretised:
            visc_residual_kernel = self.create_residual_kernel(residual_arrays, viscous_discretised, block)
            visc_residual_kernel.set_computation_name("Viscous residual")
            kernels += [visc_residual_kernel]
        #create_latex_kernel(kernels)
        # Add the kernels to the solutions
        type_of_eq.Kernels += kernels
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
        for no,array in enumerate(residual_arrays):
            expr = Eq(array, array+discretised_eq[no])
            residue_kernel.add_equation(expr)
        residue_kernel.set_grid_range(block)
        return residue_kernel
    def genral_discretisation(self, equations, block, name=None):
        """
        This discretises the central derivatives, without any special treatment to
        group the derivatives or any thing
        """
        descritised_equations = equations[:]
        cds = self.get_local_function(equations)
        if cds:
            local_kernels = {}
            if block.store_derivatives:
                for der in cds:
                    der.update_work(block)
                    ker = Kernel(block)
                    if name:
                        ker.set_computation_name("%s %s "%(name, der))
                    local_kernels[der] = ker # Reverted back
            # create a dictionary of works and kernels
            #local_kernels = {}
            work_arry_subs = {}
            for der in cds:
                self.update_range_of_constituent_relations(der, block)
                expr, local_kernels = self.traverse(der, local_kernels, block)
                expr_discretised = Eq(der.work, expr._discretise_derivative(self, block))
                work_arry_subs[expr] = der.work
                local_kernels[der].add_equation(expr_discretised)
                local_kernels[der].set_grid_range(block)
                #print "applying bcs for %s"%(expr)
                #local_kernels[der] += expr.apply_boundary_derivative_modification(block, self, der.work) # Applys the boundary kernel modifications if any
            #Apply any Boundary conditions Reverted back
            #pprint(work_arry_subs)
            for no, c in enumerate(descritised_equations):
                descritised_equations[no] = descritised_equations[no].subs(work_arry_subs)
            return local_kernels, descritised_equations
        else:
            return None, None

    def classify_equations_on_parameter(self, equations, parameter):
        from sympy import S
        containing_terms = [S.Zero for eq in equations]
        other_terms = [S.Zero  for eq in equations]
        for number, eq in enumerate(equations):
            if isinstance(eq.rhs, Add):
                for expr in eq.rhs.args:
                    if expr.has(parameter):
                        containing_terms[number] = Add(containing_terms[number],expr)
                    else:
                        other_terms[number] =Add(other_terms[number], expr)
            elif isinstance(eq.rhs, Mul):
                expr = eq.rhs
                if expr.has(parameter):
                    containing_terms[number] = Add(containing_terms[number],expr)
                else:
                    other_terms[number] =Add(other_terms[number], expr)
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
        #if CD.args[0].atoms(CentralDerivative):
        pot = postorder_traversal(CD)
        inner_cds = []
        for p in pot:
            if isinstance(p, CentralDerivative):
                inner_cds += [p]
            else:
                continue
        # Contains inner derivatives
        if len(inner_cds)>1:
            for np,cd in enumerate(inner_cds[:-1]):
                if cd.is_store and cd.work:
                    cd.is_used(True)
                    expr = expr.subs(cd, cd.work)
                    # update the kernel ranges for inner cds
                    dires = inner_cds[np+1].get_direction
                    for direction in dires:
                        kernel_dictionary[cd].set_halo_range(direction, 0, self.halotype)
                        kernel_dictionary[cd].set_halo_range(direction, 1, self.halotype)
                    #print "applying bcs for %s"%(cd)
                    #kernel_dictionary[cd] += cd.apply_boundary_derivative_modification(block, self, cd.work) # Applys the boundary kernel modifications if any (REVERTED BACK)
                elif cd.is_store:
                    # THIS raises an error when the CD(u0,x0) is not there in all derivatives ,
                    # while evaluating CD(CD(u0,x0),x1)
                    raise ValueError("NOT IMPLEMENTED THIS")
                elif not cd.is_store:# Donot store derivatives
                    raise ValueError("This dependency sgould be validated for Carpenter BC")
                    expr = expr.subs(cd, cd._discretise_derivative(self, block))
                else:
                    raise ValueError("Could not classify this")
        return expr, kernel_dictionary


from .opensbliobjects import ConstantIndexed, ConstantObject
class TemproalSolution(object):
    def __init__(self):
        self.start_kernels = []
        self.end_kernels = []
        self.kernels = []
        return
class RungeKutta(Scheme):

    """ Runge-Kutta time-stepping scheme. """

    def __init__(cls, order, constant_dt=None):
        """ Set up the Runge-Kutta stages and the coefficients.

        :arg int order: The order of accuracy of the scheme.
        """

        Scheme.__init__(cls, "RungeKutta", order)
        cls.schemetype = "Temporal"
        cls.nloops = 2
        cls.stage = Idx('stage', order)
        cls.solution_coeffs = ConstantIndexed('rkold', cls.stage)
        cls.stage_coeffs = ConstantIndexed('rknew', cls.stage)
        from .kernel import ConstantsToDeclare as CTD
        coeffs = cls.get_coefficients
        CTD.add_constant(cls.solution_coeffs, value = coeffs[cls.solution_coeffs])
        CTD.add_constant(cls.stage_coeffs, value = coeffs[cls.stage_coeffs])
        cls.solution = {}
        if constant_dt:
            raise NotImplementedError("")
        else:
            cls.constant_time_step = True
        return
    @property
    def get_coefficients(cls):
        """ Return the coefficients of the Runge-Kutta update equations.

        :returns: A dictionary of (update_equation, coefficients) pairs.
        :rtype: dict
        """

        if cls.order == 3:
            coeffs = {}
            coeffs[cls.solution_coeffs] = [Rational(1.0, 4.0), Rational(3.0, 20), Rational(3.0, 5.0)]
            coeffs[cls.stage_coeffs] = [Rational(2, 3), Rational(5, 12), Rational(3, 5)]
        return coeffs

    def __str__(cls):
        return "%s"%(cls.__class__.__name__)

    def get_local_function(cls, list_of_components):
        from .opensblifunctions import TemporalDerivative
        CD_fns = []
        for c in flatten(list_of_components):
            CD_fns += list(c.atoms(TemporalDerivative))
        return CD_fns

    def discretise(cls, type_of_eq, block):
        # We need only the equations as they contain residual residual_arrays
        if type_of_eq in cls.solution.keys():
            pass
        else:
            cls.solution[type_of_eq] = TemproalSolution()
        td_fns = cls.get_local_function(type_of_eq.equations)
        #print td_fns
        if td_fns:
            # Create a Kernel for the update ()
            old_data_sets = cls.create_old_data_sets(td_fns, block)
            new_data_sets = [eq.time_advance_array for eq in td_fns]
            zipped = zip(old_data_sets, new_data_sets )

            # create a kernel for the save equations
            kernel = cls.create_start_computations(zipped, block)
            # Add Kernel to the Solution
            cls.solution[type_of_eq].start_kernels += [kernel]
            # Create the stage and solution updates
            residuals = [eq.residual for eq in flatten(type_of_eq.equations)]
            zipped = zip(old_data_sets, new_data_sets, residuals)
            kernels = cls.create_discretisation_kernel(zipped, block)
            cls.solution[type_of_eq].kernels += kernels
            # New testing
            type_of_eq.temporalsolution = TemproalSolution()
            type_of_eq.temporalsolution.kernels += kernels
            type_of_eq.temporalsolution.start_kernels += cls.solution[type_of_eq].start_kernels
        return

    def create_discretisation_kernel(cls, zipped, block):
        solution_update_kernel = Kernel(block, "Temporal solution advancement")
        stage_update_kernel = Kernel(block, "Sub stage advancement")
        # Update the range of evaluation
        solution_update_kernel.set_grid_range(block)
        stage_update_kernel.set_grid_range(block)
        # Update the solution and stages
        if cls.constant_time_step:
            old, new = cls.constant_time_step_solution(zipped)
        solution_update_kernel.add_equation(old)
        stage_update_kernel.add_equation(new)
        solution_update_kernel.update_block_datasets(block)
        stage_update_kernel.update_block_datasets(block)
        return [stage_update_kernel, solution_update_kernel]

    def constant_time_step_solution(cls,zipped):
        dt = ConstantObject("dt")
        from .kernel import ConstantsToDeclare as CTD
        CTD.add_constant(dt)
        old = []
        new = []
        for z in zipped:
            old += [Eq(z[0], z[0] + dt*cls.solution_coeffs*z[2], evaluate=False)]
            new += [Eq(z[1], z[0] + dt*cls.stage_coeffs*z[2], evaluate=False)]
        return old, new

    def create_start_computations(cls, zipped, block):
        kernel = Kernel(block)
        kernel.computation_name = "Save equations"
        kernel.set_grid_range(block)
        for z in zipped:
            kernel.add_equation(Eq(z[0], z[1]))
        kernel.update_block_datasets(block)
        return kernel

    def create_old_data_sets(cls, equations, block):
        old_data_sets = []
        for no, eq in enumerate(flatten(equations)):
            fn = eq.time_advance_array
            old_data_sets += [block.work_array('%s_old' % fn.base.label)]
        return old_data_sets

    def generate_inner_loop(cls, kernels):
        from .algorithm import DoLoop
        rkloop = DoLoop(cls.stage)
        rkloop.add_components(kernels)
        return rkloop
