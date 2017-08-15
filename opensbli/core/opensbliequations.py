from .opensbliobjects import DataSet, ConstantObject
from .opensblifunctions import *
from .kernel import Kernel
from sympy import flatten, preorder_traversal
from sympy import Equality
import copy
from sympy.tensor.array import MutableDenseNDimArray
class Discretisation(object):
    """
    This should contain the following
    1. Required_datasets these are the required data sets used in the simulations
    2. Required_constants these are the required constants for the SimulationEquations
    3. Local_evaluations this returns a list of evluations that should be performed for
    the Discretisation of any equations
    """
    @property
    def required_datasets(cls):
        """By the time this function is called all the functions such as
        KD, LC, DOT etc.. should be evaluated
        """
        objs = []
        for eq in flatten(cls.equations):
            objs += list(eq.atoms(DataSet))
        objs = set(objs)
        return objs
    @property
    def required_constants(cls):
        constants = []
        for eq in flatten(cls.equations):
            constants += list(eq.atoms(ConstantObject))
        constants = set(constants)
        return constants
    @property
    def required_functions_local(cls):
        cds = []
        for eq in flatten(cls.equations):
            pot = preorder_traversal(eq)
            for p in pot:
                if isinstance(p, Function):
                    cds +=[p]
                    pot.skip()
                else:
                    continue
        fns = set(cds)
        return fns
    @property
    def required_functions(cls):
        """
        Returns all the functions to be evaluated
        """
        fns = cls.required_functions_local
        allfns = flatten([[fn,fn.required_functions] for fn in fns])
        allfns = set(allfns)
        return allfns
    def apply_CentralDerivative(cls):
        """
        This creates a list of central derivatives with
        """
        fns = cls.required_functions
        for fn in fns:
            pass
        return
    def _sanitise_equations(cls, equation):
        fns = []
        replacements ={}
        if isinstance(equation, list):
            for no,eq in enumerate(equation):
                fns += list(eq.atoms(Function))
        else:
            fns += list(equation.atoms(Function))
        for fn in fns:
            if not fn.is_homogeneous:
                replacements[fn] = fn._sanitise
        # modify the original equations
        if isinstance(equation, list):
            for no,eq in enumerate(equation):
                equation[no] = equation[no].xreplace(replacements)
        else:
            equation = equation.xreplace(replacements)
        return equation

    def apply_metrics(cls, metriceqclass):
        out = []
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

class OpenSBLIEquation(Equality):
    is_Equality = True
    def __new__(cls, lhs, rhs):
        ret = super(OpenSBLIEquation, cls).__new__(cls, lhs, rhs)
        ret.is_vector = False
        return ret
    def _eval_evalf(self, prec):
        return self.func(*[s._evalf(prec) for s in self.args])
    def set_vector(cls, component_number):
        cls.is_vector = True
        cls.vector_component = component_number
        return
    @property
    def _sanitise_equation(cls):
        fns = []
        replacements ={}
        fns += list(cls.atoms(Function))

        for fn in fns:
            if not fn.is_homogeneous:
                replacements[fn] = fn._sanitise
        eq = cls
        return eq.xreplace(replacements)

class Solution(object):
    def __init__(self):
        # Kernels would be spatial kernels
        self.Kernels = []
        self.constituent_relations_kernels = {}
        # Kernels for boundary
        self.boundary_kernels = []
        return
class OpenSBLIExpression(Expr):
    """
    This represents each and every Discretisation expression
    for example, to evaluate CD((p+rhoE)*u0, x0)
    will have two expressions
    a. to evaluate work = (p+rhoE)*u0
    b. to evaluate next_work = CD(work, x0)

    """
    def __new__(cls, expr):
        ret = Expr.__new__(cls, expr)
        return ret
    @property
    def as_expr(cls):
        return cls.args[0]

class SimulationEquations(Discretisation, Solution):
    """
    To proceed
    We will create solution equations, these are nothing but the
    steps required to solve the equations. For example,
    to solve say continuity equation, we already have (rhou0:ndim)
    so the equations would be
    Eq(CD(rhou0,x0), ITS descritised formula)
    similarly, for the momentum equation
    or, the best way would be create an object for evaluating for each and every
    term of the equation CD / WD, TD
    CD(rhou0,x0,x0) --> is an evaluation object already have function (CD).
    It should give you requires

    1. mark all the derivatives to be stored or to be evaluated
    on the fly.
    2. If stored and contain any dependant derivatives then these should be replaced
    by the inner derivatives
        example, CD(u0,x0,x1) --> CD(CD(u0,x0),x1)
        CD(u0,x0,x0) should not be replaced

    This should have the following functions
    a. converting the DataObjects to DataSet
    b. Applying the functions (KD, LC etc..)
    c. Applying derivatives (value)
    a. required functions (spatial and temporal)

    """
    def __new__(cls, order = None, **kwargs):
        ret = super(SimulationEquations,cls).__new__(cls)
        if order:
            ret.order = order
        else:
            ret.order = 0
        ret.equations = []
        ret.kwargs = kwargs
        return ret

    def add_equations(cls, equation):
        #equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            local = []
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq.lhs, eq.rhs)
                eq.set_vector(no)
                local += [eq]
            cls.equations += [local]
        else:
            equation = OpenSBLIEquation(equation.lhs, equation.rhs)
            cls.equations += [equation]
        return

    def add_constituent_relations(cls, constituent_relations):
        """
        Adds the constituent relations we will make a deep copy so that the equations
        can be reused else where
        """
        #cls.constituent_relations = constituent_relations
        cls.constituent_relations_dictionary = cls.convert_to_dictionary(constituent_relations.equations)
        return

    def convert_to_dictionary(cls, relations):
        output_dictionary = {}
        for r in relations:
            output_dictionary[r.lhs] =  r.rhs
        return output_dictionary

    def create_residual_arrays(cls, block):
        for no, eq in enumerate(flatten(cls.equations)):
            if not hasattr(eq, 'residual'):
                eq.residual = DataSetBase('Residual%d'%no)[[0 for i in range(block.ndim)]]
        return
    def zero_residuals_kernel(cls, block):
        """Evaluate all the residuals to zero
        """
        kernel = Kernel(block, computation_name= "Zeroing residuals")
        eqs = [Eq(eq.residual, 0) for eq in flatten(cls.equations)]
        kernel.add_equation(eqs)
        kernel.set_grid_range(block)
        cls.Kernels += [kernel]
        return


    @property
    def get_required_constituents(self):
        arrays = []
        for eq in flatten(self.equations):
            arrays += list(eq.atoms(DataSet))
        arrays = set(arrays) 
        return arrays

    def spatial_discretisation(cls, schemes, block):
        """
        Discretisation will be done by
        a. Update the work arrays (based on location) for the CD, WD, TD --> irrespective of Type
        b. Set the evaluation range of each of the workbases to be that of the grid, irrespective of CD, WD or TD
        c. Descritise each and every CD, WD, TD
        """
        # Instantiate the solution class
        (Solution,cls).__init__(cls)
        # Create the residual array for the equations
        cls.create_residual_arrays(block)
        cls.zero_residuals_kernel(block)

        #cls.descritsed_equations = copy.copy(cls.equations)
        spatialschemes = []
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation
        cls.constituent_evaluations = {}
        crs = block.get_constituent_equation_class
        cr_dictionary = {}
        for cr in crs:
            cr_dictionary.update(cr.get_relations_dictionary)
        cls.requires = {}
        for no,sc in enumerate(spatialschemes):
            cls.constituent_evaluations[sc] = schemes[sc].discretise(cls, block)
            for key, value in cls.constituent_evaluations[sc].iteritems():
                if key in cr_dictionary.keys():
                    if key in cls.constituent_relations_kernels:
                        cls.constituent_relations_kernels[key].merge_halo_range(value.halo_ranges)
                    else:
                        cls.constituent_relations_kernels[key] = value
                        cls.constituent_relations_kernels[key].add_equation(cr_dictionary[key])
                else:
                    #raise ValueError("Constituent relation is not found for %s"%key)
                    cls.requires[key] = value
        missing_CR_datasets = cls.get_required_constituents.difference(cls.constituent_relations_kernels.keys())
        for dset in missing_CR_datasets:
            # Evaluation of missing dataset is required
            if dset in cr_dictionary.keys():
                for kernel in cr_dictionary[dset].kernels:
                    cls.constituent_relations_kernels[kernel.equations[0].lhs] = kernel
        cls.process_kernels(block)

        return
    def process_kernels(cls, block):
        for key,kernel in cls.constituent_relations_kernels.iteritems():
            kernel.update_block_datasets(block)
        for kernel in cls.Kernels:
            kernel.update_block_datasets(block)
        #TODO similarly  update the rational constants and constants
        """
        Constants should have the attributes
        a. is_input = True or False
        b. dtype = int or double or float
        c. if not input the value
        """
        return
    @property
    def sort_constituents(cls):
        """
        Sort the constituent relations thinking requires is Known
        """
        input_order = []
        for a in cls.requires.keys():
            if isinstance(a, DataSet):
                input_order += [a.base]
            else:
                input_order +=  [a]

        dictionary = {}
        for key, value in cls.constituent_relations_kernels.iteritems():
            dictionary[key.base] = value
        order_of_evaluation = cls.sort_dictionary(input_order, dictionary)
        ordered_kernels = []
        for o in order_of_evaluation:
            ordered_kernels += [dictionary[o]]
        return ordered_kernels
    def sort_dictionary(cls, order, new_dictionary):
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
        order = flatten(order + [a.base for a in flatten(cls.time_advance_arrays)])
        order = list(set(order))
        # store the length of order
        input_order = len(order)
        key_list = [key for key in dictionary.keys() if key not in order]
        requires_list = ([dictionary[key].required_data_sets for key in key_list])

        zipped = zip(key_list, requires_list)
        # Breaks after 1000 iterations
        iter_count = 0
        while key_list:
            iter_count = iter_count+1
            order += [x for (x, y) in zipped if all(req in order for req in y)]
            key_list = [key for key in dictionary.keys() if key not in order]
            requires_list = [dictionary[key].required_data_sets for key in key_list]
            zipped = zip(key_list, requires_list)
            if iter_count > 1000:
                print("Exiting because i cannot classify the following")
                print("Already sorted are")
                pprint(order)
                pprint([srepr(o) for o in order])
                #print("Trying to sort the required for")
                #pprint(evaluations[key].lhs)
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
                    local += [td.time_advance_array  for td in c1.atoms(TemporalDerivative)]
                TD_fns += [local]
            else:
                TD_fns += [td.time_advance_array  for td in c.atoms(TemporalDerivative)]
        return TD_fns
    @property
    def algorithm_location(cls):
        return True
    @property
    def all_spatial_kernels(cls):
        return cls.sort_constituents + cls.Kernels

class ConstituentRelations(Discretisation, Solution):
    def __new__(cls):
        ret = super(ConstituentRelations,cls).__new__(cls)
        ret.equations = []
        ret.vector_number = 0 # Used later a place holder for multiple vectors
        ret.order = None # This means there is no location for this explicitly in the algorithm
        return ret

    def add_equations(cls, equation):
        #equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq.lhs, eq.rhs)
                #eq.set_vector(no)
                cls.equations += [eq]
        else:
            equation = OpenSBLIEquation(equation.lhs, equation.rhs)
            cls.equations += [equation]
        return

    def create_residual_arrays(cls):
        for eq in flatten(cls.equations):
            eq.residual = eq.lhs
        return

    def spatial_discretisation(cls, schemes, block):
        """
        Discretisation will be done by
        a. Update the work arrays (based on location) for the CD, WD, TD --> irrespective of Type
        b. Set the evaluation range of each of the workbases to be that of the grid, irrespective of CD, WD or TD
        c. Descritise each and every CD, WD, TD
        """
        # Instantiate the solution class
        (Solution,cls).__init__(cls)
        # Create the residual array for the equations
        cls.create_residual_arrays()

        #cls.descritsed_equations = copy.copy(cls.equations)
        spatialschemes = []
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        # Perform spatial Discretisation if any in constituent relations evaluation
        cls.constituent_evaluations = {}
        equations = cls.equations

        for eq in equations:
            cls.equations = [eq]
            for sc in spatialschemes:
                cls.constituent_evaluations[sc] = schemes[sc].discretise(cls, block)
            eq.kernels = cls.Kernels[:]

            cls.Kernels = []
        cls.equations = equations
        for sc in spatialschemes:
            cls.constituent_evaluations[sc] = schemes[sc].discretise(cls, block)
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

class NonSimulationEquations():
    """ Dummy place holder for all the equations that are not simulated but needs to be evaluated
    e.g, metrics or diagnostics or Statistics, """
    pass


class DiagnosticsEquations(NonSimulationEquations):

    def __new__(cls):
        # The order for this would be >100
        return
