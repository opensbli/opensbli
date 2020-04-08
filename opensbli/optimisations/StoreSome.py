from sympy import flatten, simplify, symbols
from opensbli.core.opensbliobjects import ConstantObject, CoordinateObject, DataObject
from opensbli.core.grid import GridVariable
from opensbli.core.opensblifunctions import CentralDerivative
from opensbli.core.kernel import Kernel
from opensbli.schemes.spatial import Central
from opensbli.equation_types.opensbliequations import OpenSBLIEq, SimulationEquations


class StoreSome(Central):
    """ Low-storage algorithms to reduce memory intensity and the number of global storage arrays.
        S.P. Jammy et al. Journal of Computational Science. Vol 36, September 2019 10.015."""

    def __init__(self, order, der_fns_to_store, level=1):
        """ Set up the scheme.
        :arg int order: The order of accuracy of the scheme."""
        Central.__init__(self, order)
        self.fns = der_fns_to_store
        return

    def derivatives_to_store(self, coordinates, block):
        """ Creates CentralDerivative objects of the derivatives to be stored."""
        data_objects = flatten([symbols('%s' % self.fns, **{'cls': DataObject})])
        data_sets = [block.location_dataset(str(d)) for d in data_objects]
        coords = [c for c in coordinates if not c.get_coordinate_type()]
        self.derivatives_to_store = []
        for a in data_sets:
            for b in coords:
                self.derivatives_to_store += [CentralDerivative(a, b)]
        return

    def discretise(self, type_of_eq, block):
        """ This is the main calling function from opensbli equations.spatial_discretisation, which is called from block.discretise."""
        self.set_halos(block)
        if isinstance(type_of_eq, SimulationEquations):
            # Simulation equations are always solved as sbli_rhs_discretisation currently
            self.sbli_rhs_discretisation(type_of_eq, block)
            return self.required_constituent_relations
        else:
            discretised_eq = self.SS(type_of_eq.equations, block, None, group=False)
            if discretised_eq:
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

    def sbli_rhs_discretisation(self, type_of_eq, block):
        """ Function that performs the discretisation to the convective and viscous terms."""
        block.reset_work_index
        equations = flatten(type_of_eq.equations)
        residual_arrays = [eq.residual for eq in equations]
        equations = [e._sanitise_equation for e in equations]
        coordinates = set()
        for e in equations:
            coordinates = coordinates.union(e.atoms(CoordinateObject))
        self.derivatives_to_store(coordinates, block)
        self.local_kernels = {}
        self.required_constituent_relations = {}
        subs_dict = {}
        for no, der in enumerate(self.derivatives_to_store):
            der.update_work(block)
            subs_dict[der] = der.work
            ker = Kernel(block)
            ker.set_computation_name("Derivative evaluation %s " % (der))
            self.update_range_of_constituent_relations(der, block)
            v = der
            expr = OpenSBLIEq(v.work, simplify(v._discretise_derivative(self, block)))
            # pprint(expr)
            ker.add_equation(expr)
            ker.set_grid_range(block)
            self.local_kernels[v] = ker
            self.derivatives_to_store[no] = v

        # Get all the derivatives and traverse them to update the range of
        # already evaluated derivatives
        for d in self.get_local_function(equations):
            if d.atoms(CentralDerivative).intersection(self.derivatives_to_store):
                expr, self.local_kernels = self.traverse(d, self.local_kernels, block)
        equations = [e.subs(subs_dict) for e in equations]
        classify_parameter = ConstantObject("Re")
        rhs_eq = [e.rhs for e in equations]
        discretised_eq = [OpenSBLIEq(x, y) for x, y in zip(residual_arrays, rhs_eq)]

        viscous, convective = self.classify_equations_on_parameter(discretised_eq, classify_parameter)
        # Apply to the convective terms
        convective = [OpenSBLIEq(x, y) for x, y in zip(residual_arrays, convective)]
        convective_equations = self.SS(convective, block, 'Convective')
        # Apply to the viscous terms
        viscous = [OpenSBLIEq(x, x+y) for x, y in zip(residual_arrays, viscous)]
        viscous_equations = self.SS(viscous, block, 'Viscous')
        if convective_equations or viscous_equations:
            for ker in self.local_kernels:
                eval_ker = self.local_kernels[ker]
                type_of_eq.Kernels += [eval_ker]
        if convective_equations:
            convective_kernel = Kernel(block, computation_name="Convective terms")
            convective_kernel.set_grid_range(block)
            for eq in convective_equations:
                convective_kernel.add_equation(eq)
            convective_kernel.update_block_datasets(block)
            type_of_eq.Kernels += [convective_kernel]
        if viscous_equations:
            viscous_kernel = Kernel(block, computation_name="Viscous terms")
            viscous_kernel.set_grid_range(block)
            for eq in viscous_equations:
                viscous_kernel.add_equation(eq)
            viscous_kernel.update_block_datasets(block)
            type_of_eq.Kernels += [viscous_kernel]
        block.reset_work_index
        return self.required_constituent_relations

    def SS(self, equations, block, equation_type, group=True, level=1):
        """ Generates the GridVariables to store the local derivatives. Also sorts the
        derivatives to re-access the same arrays consecutively."""
        discrete_equations = flatten(equations)[:]
        cds = self.get_local_function(flatten(equations))
        grid_variable_evaluations = []
        if cds:
            gridvars = [GridVariable('localeval_%d' % i) for i in range(len(cds))]
            # Sort by grouping variables
            if group:
                if equation_type is 'Convective':
                    cds = sorted(cds, key=lambda x: x.args[1].direction)
                elif equation_type is 'Viscous':
                    cds = sorted(cds, key=lambda x: str(x.args[0]))
            for der in cds:
                self.update_range_of_constituent_relations(der, block)
                if level == 1:
                    gv = gridvars.pop(0)
                    grid_variable_evaluations += [OpenSBLIEq(gv, simplify(der._discretise_derivative(self, block)))]
                    for no, c in enumerate(discrete_equations):
                        discrete_equations[no] = discrete_equations[no].subs(der, gv)
            return grid_variable_evaluations+discrete_equations
        else:
            return None
