"""@brief Bsae of discretisaion of equations
   @author Satya Pramod Jammy
   @contributors David Lusher
   @details base classes for different type of equations used in opensbli
"""

from opensbli.core.opensbliobjects import DataSet, ConstantObject, DataSetBase, DataObject
from opensbli.core.opensblifunctions import TemporalDerivative
from sympy import flatten, preorder_traversal
from sympy import Equality, Function, pprint, srepr


class Discretisation(object):
    """Contains the functions used in various equation classes of OpenSBLI, to perform
    discretisation. """

    @property
    def required_datasets(cls):
        """Returns the set of all datasets in the equations of a class.

        :return: a set of DataSet
        :rtype: set(DataSet) """
        objs = []
        for eq in flatten(cls.equations):
            objs += list(eq.atoms(DataSet))
        objs = set(objs)
        return objs

    @property
    def required_constants(cls):
        """Returns the set of all constants in the equations.

        :return: a set of constants
        :rtype: set(ConstantObject) """
        constants = []
        for eq in flatten(cls.equations):
            constants += list(eq.atoms(ConstantObject))
        constants = set(constants)
        return constants

    @property
    def required_functions_local(cls):
        """Helper funciton to evaluate the required finctions."""
        cds = []
        for eq in flatten(cls.equations):
            pot = preorder_traversal(eq)
            for p in pot:
                if isinstance(p, Function):
                    cds += [p]
                    pot.skip()
                else:
                    continue
        fns = set(cds)
        return fns

    @property
    def required_functions(cls):
        """Returns all the functions in the equations.

        :return: a set of functions
        :rtype: set(CentralDerivative, WenoDerivative, TemporalDerivative, etc.) """
        fns = cls.required_functions_local
        allfns = flatten([[fn, fn.required_functions] for fn in fns])
        allfns = set(allfns)
        return allfns

    def _sanitise_equations(cls, equation):
        """Sanitises the equation. Finds the non homogenous functions and replace them with
        a funciton(function) approach.

        .. note:
            **Example**\n
            CentralDerivative(u0, x0, x1) is transformed to CetralDerivative(CentralDerivative(u0, x0), x1)
        """
        fns = []
        replacements = {}
        if isinstance(equation, list):
            for no, eq in enumerate(equation):
                fns += list(eq.atoms(Function))
        else:
            fns += list(equation.atoms(Function))
        for fn in fns:
            if not fn.is_homogeneous:
                replacements[fn] = fn._sanitise
        # modify the original equations
        if isinstance(equation, list):
            for no, eq in enumerate(equation):
                equation[no] = equation[no].xreplace(replacements)
        else:
            equation = equation.xreplace(replacements)
        return equation

    def apply_metrics(cls, metriceqclass):
        """Applies the metric transformation for the equations in the class using the metric class provided

        :param MetricsEquation metriceqclass: see :class:`.MetricsEquation`
        :return: None"""
        out = []
        if len(flatten(cls.equations)) == 0:
            raise ValueError("%s class does not have equations to apply metrics to." % cls.__class__.__name__)
        for eq in cls.equations:
            if isinstance(eq, list):
                out_inner = []
                for eq1 in eq:
                    out_inner += [metriceqclass.apply_transformation(eq1)]
                out += [out_inner]
            else:
                out += [metriceqclass.apply_transformation(eq)]
        cls.equations = out
        return

    def write_latex(self, latex):
        """Writes the latex form of the equations in the class to the file

        :param LatexWriter latex"""
        latex.write_string('The euqations are %s' % type(self).__name__)
        for eq in flatten(self.equations):
            latex.write_expression(eq)
        return

    def convert_to_dictionary(self, relations):
        """Helper function to convert a list of relations to a dictionary"""
        output_dictionary = {}
        for r in relations:
            output_dictionary[r.lhs] = r.rhs
        return output_dictionary

    def add_equations(self, equation):
        """Adds the equations to the class

        :param equation: a list of equations or a single equation to be added to the class"""
        if isinstance(equation, list):
            local = []
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq.lhs, eq.rhs)
                eq.set_vector(no)
                local += [eq]
            self.equations += [local]
        else:
            equation = OpenSBLIEquation(equation.lhs, equation.rhs)
            self.equations += [equation]
        return

    @property
    def evaluated_datasets(cls):
        """This should be over written in different classes, provided here to give a default option"""
        return set()


class OpenSBLIEquation(Equality):
    """A wrapper around SymPy's Equality class to provide more flexibility in the future"""
    is_Equality = True

    def __new__(cls, lhs, rhs, **options):
        ret = super(OpenSBLIEquation, cls).__new__(cls, lhs, rhs, **options)
        return ret

    def set_vector(cls, component_number):
        cls.is_vector = True
        cls.vector_component = component_number
        return

    def convert_to_datasets(self, block):
        replacements = {}
        for d in self.atoms(DataObject):
            replacements[d] = block.location_dataset(d)
        args = [a.subs(replacements) for a in self.args]
        return type(self)(*args)

    @property
    def required_datasets(self):
        """Not required"""
        return self.lhs_datasetbases.union(self.rhs_datasetbases)

    @property
    def lhs_datasetbases(self):
        return self.lhs.atoms(DataSetBase)

    @property
    def rhs_datasetbases(self):
        return self.rhs.atoms(DataSetBase)

    @property
    def _sanitise_equation(cls):
        fns = []
        replacements = {}
        fns += list(cls.atoms(Function))

        for fn in fns:
            if not fn.is_homogeneous:
                replacements[fn] = fn._sanitise
        eq = cls
        return eq.xreplace(replacements)


OpenSBLIEq = OpenSBLIEquation


class Solution(object):
    """An object to store the symbolic kernels generated while applying the numerical method to the equations"""

    def __init__(self):
        # Kernels would be spatial kernels
        self.Kernels = []
        self.constituent_relations_kernels = {}
        # Kernels for boundary
        self.boundary_kernels = []
        return


class SimulationEquations(Discretisation, Solution):
    """Class for the simulation equations. This performs the discretisation of the equations.

    :param int order: priority in the algorithm if multiple simulation equations exitst
    """

    def __new__(cls, order=None, **kwargs):
        ret = super(SimulationEquations, cls).__new__(cls)
        if order:
            ret.order = order
        else:
            ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        return ret

    def add_constituent_relations(cls, constituent_relations):
        """Convert the constituent relations passed to a dictionary for easier access in discretisation

        :param ConstituentRelations constituent_relations: the constituent relations class on the block"""

        cls.constituent_relations_dictionary = cls.convert_to_dictionary(constituent_relations.equations)
        return

    def create_residual_arrays(cls, block):
        """Creates the residual datasets for each of the simulation equations.

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        for no, eq in enumerate(flatten(cls.equations)):
            if not hasattr(eq, 'residual'):
                eq.residual = block.location_dataset('Residual%d' % no)
        return

    # def zero_residuals_kernel(cls, block):
    #     """A symbolic kernel to equate all the residual arrays of the equations to zero

    #     :param SimulationBlock block: the block on which the equations are solved
    #     :return: None """
    #     from opensbli.core.kernel import Kernel
    #     kernel = Kernel(block, computation_name="Zeroing residuals")
    #     eqs = [OpenSBLIEq(eq.residual, 0) for eq in flatten(cls.equations)]
    #     kernel.add_equation(eqs)
    #     kernel.set_grid_range(block)
    #     cls.Kernels += [kernel]
    #     return

    @property
    def get_required_constituents(self):
        arrays = []
        for eq in flatten(self.equations):
            arrays += list(eq.atoms(DataSet))
        arrays = set(arrays)
        return arrays

    def spatial_discretisation(cls, block):
        """Applies the spatial discretisation of the equations by calling the discretisation of each spatial 
        scheme provided on the block

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """

        # Instantiate the solution class
        (Solution, cls).__init__(cls)
        # Create the residual array for the equations
        cls.create_residual_arrays(block)
        # Kernel to make the residuals zero
        # cls.zero_residuals_kernel(block)

        spatialschemes = []
        # Get the schemes on the block
        schemes = block.discretisation_schemes
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation
        cls.constituent_evaluations = {}
        crs = block.get_constituent_equation_class
        missing_CR_datasets = []
        cr_dictionary = {}
        for cr in crs:
            cr_dictionary.update(cr.get_relations_dictionary)
        cls.requires = {}
        for no, sc in enumerate(spatialschemes):
            cls.constituent_evaluations[sc] = schemes[sc].discretise(cls, block)
            for key, value in cls.constituent_evaluations[sc].iteritems():
                if key in cr_dictionary.keys():
                    if key in cls.constituent_relations_kernels:
                        if cr_dictionary[key].kernels:
                            raise NotImplementedError('')
                        cls.constituent_relations_kernels[key].merge_halo_range(value.halo_ranges)
                    else:
                        cls.constituent_relations_kernels[key] = value
                        if cr_dictionary[key].kernels:
                            missing_CR_datasets += [key]
                        else:
                            cls.constituent_relations_kernels[key].add_equation(cr_dictionary[key])
                else:
                    # raise ValueError("Constituent relation is not found for %s"%key)
                    cls.requires[key] = value
        missing_CR_datasets += list(cls.get_required_constituents.difference(cls.constituent_relations_kernels.keys()))
        for dset in missing_CR_datasets:
            # Evaluation of missing dataset is required
            if dset in cr_dictionary.keys():
                for kernel in cr_dictionary[dset].kernels:
                    cls.constituent_relations_kernels[kernel.equations[-1].lhs] = kernel
        cls.process_kernels(block)
        return

    def process_kernels(cls, block):
        """A function to update some dependant parameters of each kernel

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        for key, kernel in cls.constituent_relations_kernels.iteritems():
            kernel.update_block_datasets(block)
        for kernel in cls.Kernels:
            kernel.update_block_datasets(block)
        return

    def sort_constituents(cls, block):
        """Sort the constituent relation kernels
        """
        input_order = []
        for a in cls.requires.keys():
            if isinstance(a, DataSet):
                input_order += [a.base]
            else:
                input_order += [a]

        dictionary = {}
        for key, value in cls.constituent_relations_kernels.iteritems():
            dictionary[key.base] = value
        order_of_evaluation = cls.sort_dictionary(input_order, dictionary, block)
        ordered_kernels = []
        for o in order_of_evaluation:
            ordered_kernels += [dictionary[o]]
        return ordered_kernels

    @property
    def evaluated_datasets(cls):
        return set([c.base for c in flatten(cls.time_advance_arrays)])

    def sort_dictionary(cls, order, new_dictionary, block):
        """ Sort the evaluations based on the requirements of each term. For example, if we have
        the primitive variables p, u0, u1, and T, then the pressure p may depend on the velocity u0 and u1, and T may depend on p,
        so we need this be evaluate in the following order: u0, u1, p, T.


        :arg list order: The list of already sorted terms or Known terms
        :arg evaluations: The evaluation information, containing dependency information.
        :arg typef: The type of term to sort.
        :returns: A list of ordered terms.
        :rtype: list
        """
        dictionary = new_dictionary
        # reverse_dictionary = {}
        order = flatten(order + list(block.known_datasets))
        order = list(set(order))
        # store the length of order
        input_order = len(order)
        key_list = [key for key in dictionary.keys() if key not in order]
        requires_list = ([dictionary[key].rhs_datasetbases for key in key_list])
        zipped = zip(key_list, requires_list)
        # Breaks after 1000 iterations
        iter_count = 0
        while key_list:
            iter_count = iter_count+1
            order += [x for (x, y) in zipped if all(req in order for req in y)]
            key_list = [key for key in dictionary.keys() if key not in order]
            requires_list = [dictionary[key].rhs_datasetbases for key in key_list]
            zipped = zip(key_list, requires_list)
            if iter_count > 1000:
                print("Exiting because i cannot classify the following")
                print("Already sorted are")
                pprint(order)
                pprint([srepr(o) for o in order])
                # print("Trying to sort the required for")
                # pprint(evaluations[key].lhs)
                print("It requires")
                pprint([req for req in requires_list[0]])
                print("Sorted")
                pprint([(req, req in order) for req in requires_list[0]])
                raise ValueError("Exiting sort evaluations ")
        order = order[input_order:]
        return order

    def temporal_discretisation(cls, schemes, block):
        """
        This should return a temporal solution class
        """

        return

    def apply_boundary_conditions(cls, block):
        arrays = cls.time_advance_arrays
        kernels = block.apply_boundary_conditions(arrays)
        cls.boundary_kernels += kernels
        return

    @property
    def time_advance_arrays(cls):
        TD_fns = []
        for c in cls.equations:
            if isinstance(c, list):
                local = []
                for c1 in c:
                    local += [td.time_advance_array for td in c1.atoms(TemporalDerivative)]
                TD_fns += [local]
            else:
                TD_fns += [td.time_advance_array for td in c.atoms(TemporalDerivative)]
        return TD_fns

    @property
    def algorithm_location(cls):
        return True

    def all_spatial_kernels(cls, block):
        return cls.sort_constituents(block) + cls.Kernels


class ConstituentRelations(Discretisation, Solution):
    """Class for the ConstituentRelations to performs the discretisation of the equations"""
    def __new__(cls):
        ret = super(ConstituentRelations, cls).__new__(cls)
        ret.equations = []
        ret.vector_number = 0  # Used later a place holder for multiple vectors
        ret.order = None  # This means there is no location for this explicitly in the algorithm
        return ret

    def create_residual_arrays(cls):
        """The Residual arrays are currently the lhs of the constituent relations"""
        for eq in flatten(cls.equations):
            eq.residual = eq.lhs
        return

    def spatial_discretisation(cls, block):
        """Applies the spatial discretisation of the equations by calling the discretisation of each spatial
        scheme provided on the block

        :param SimulationBlock block: the block on which the equations are solved
        :return: None """
        # Instantiate the solution class
        (Solution, cls).__init__(cls)
        # Create the residual array for the equations
        cls.create_residual_arrays()

        spatialschemes = []
        # Get the schemes on the block
        schemes = block.discretisation_schemes
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation if any in constituent relations evaluation
        cls.constituent_evaluations = {}
        equations = cls.equations

        for eq in flatten(equations):
            cls.equations = [eq]
            for sc in spatialschemes:
                cls.constituent_evaluations[sc] = schemes[sc].discretise(cls, block)
            eq.kernels = cls.Kernels[:]

            cls.Kernels = []
        cls.equations = equations
        return

    @property
    def get_relations_dictionary(cls):
        relations_dictionary = {}
        for eq in flatten(cls.equations):
            relations_dictionary[eq.lhs] = eq
        return relations_dictionary

    def apply_boundary_conditions(cls, block):
        pass
        return


class NonSimulationEquations(Discretisation):
    """ Dummy place holder for all the equations that are not simulated but needs to be evaluated
    e.g, metrics or post process equations or user defined kernels"""
    pass
