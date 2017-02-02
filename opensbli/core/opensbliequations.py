from .opensbliobjects import DataSet, ConstantObject
from .opensblifunctions import *
from sympy import flatten, preorder_traversal
from sympy import Equality
import copy
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


class OpenSBLIEquation(Equality):
    def __new__(cls, equation):
        ret = super(OpenSBLIEquation, cls).__new__(cls, equation.lhs, equation.rhs)
        ret.is_vector = False
        return ret
    def set_vector(cls, component_number):
        cls.is_vector = True
        cls.vector_component = component_number
        return

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
        equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            local = []
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq)
                eq.set_vector(no)
                local += [eq]
            cls.equations += [local]
        else:
            equation = OpenSBLIEquation(equation)
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

    def create_residual_arrays(cls):
        for no, eq in enumerate(flatten(cls.equations)):
            if not hasattr(eq, 'residual'):
                eq.residual = DataSet('Residual%d'%no)
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
                    cls.requires[key] = value
        cls.process_kernels(block)
        return
    def process_kernels(cls, block):
        pprint(cls.constituent_relations_kernels)
        for key,kernel in cls.constituent_relations_kernels.iteritems():
            kernel.update_block_datasets(block)
        for kernel in cls.Kernels:
            kernel.update_block_datasets(block)
        return
    @property
    def sort_constituents(cls):
        """
        Sort the constituent relations thinking requires is Known
        """
        order_of_evaluation = cls.sort_dictionary(list(cls.requires.keys()), cls.constituent_relations_kernels)
        ordered_kernels = []
        for o in order_of_evaluation:
            ordered_kernels += [cls.constituent_relations_kernels[o]]
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
        reverse_dictionary = {}
        order = flatten(order + cls.time_advance_arrays)
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
            #requires_list = flatten(requires_list)
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
        equation = cls._sanitise_equations(equation)
        if isinstance(equation, list):
            for no, eq in enumerate(equation):
                eq = OpenSBLIEquation(eq)
                eq.set_vector(no)
                cls.equations += [eq]
        else:
            equation = OpenSBLIEquation(equation)
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

class NonSimulationEquations(Discretisation, Solution):
    """ Dummy place holder for all the equations that are not simulated but needs to be evaluated
    e.g, metrics or diagnostics or Statistics, """
    pass


class DiagnosticsEquations(NonSimulationEquations):

    def __new__(cls):
        # The order for this would be >100
        return
from .parsing import Equation
class MetricsEquation(Discretisation):

    def __new__(cls, ndim, coordinate_symbol, parameters, max_order):
        ret = super(MetricsEquation,cls).__new__(cls)
        ret.ndim = ndim
        ret.stretching_metric = [param[0] for param in parameters] # Get whether true or false fr streching
        ret.curvilinear_metric = [param[1] for param in parameters] # Get T/F for curvilinear
        cart = CoordinateObject('%s_i'%(coordinate_symbol))
        curv = CoordinateObject('%s_i'%('xi'))
        cartesian_coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        curvilinear_coordinates = [curv.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        ret.curvilinear_coordinates = curvilinear_coordinates
        ret.cartesian_coordinates = cartesian_coordinates
        ret.cart_to_curvilinear_functions = Matrix(ndim, ndim,lambda i,j:curvilinear_coordinates)
        ret.curvilinear_to_cart_functions = Matrix(ndim, ndim,lambda i,j:cartesian_coordinates)
        ret.update_curvilinear_function(ret.curvilinear_metric)
        ret.update_stretching_function(ret.stretching_metric)
        ret.update_functions()
        pprint(ret.curvilinear_metric)
        pprint(ret.stretching_metric)

        #pprint(ret.cart_to_curvilinear_functions)
        #pprint(ret.curvilinear_to_cart_functions)
        #ret.generate_FD(coordinate_symbol)
        ret.generate_jacobians(coordinate_symbol, max_order)
        ret.order = -1 # This is before the time loop in the algorithm
        return ret

    def update_functions(cls):
        for d in range(cls.ndim):
            for d1 in range(cls.ndim):
                list1 = cls.cart_to_curvilinear_functions[d,d1]
                if not any(list1):
                    if d==d1:
                        #cls.cart_to_curvilinear_functions[d,d1] = cls.curvilinear_coordinates[d]
                        cls.cart_to_curvilinear_functions[d,d1] = 1
                    else:
                        cls.cart_to_curvilinear_functions[d,d1] = 0
                else:
                    inds = [ind for ind,val in enumerate(list1) if not val]
                    list1 = [v for i, v in enumerate(list1) if i not in inds]
                    #cls.cart_to_curvilinear_functions[d,d1] = cls.cartesian_coordinates[d](*list1)
                    cls.cart_to_curvilinear_functions[d,d1] = 1
                #  curvilinear to cartesian
                list1 = cls.curvilinear_to_cart_functions[d,d1]
                if not any(list1):
                    if d==d1:
                        #cls.curvilinear_to_cart_functions[d,d1] = cls.cartesian_coordinates[d]
                        cls.curvilinear_to_cart_functions[d,d1] = 1
                    else:
                        cls.curvilinear_to_cart_functions[d,d1] = 0
                else:
                    inds = [ind for ind,val in enumerate(list1) if not val]
                    list1 = [v for i, v in enumerate(list1) if i not in inds]
                    #cls.curvilinear_to_cart_functions[d,d1] = cls.curvilinear_coordinates[d](*list1)
                    #cls.curvilinear_to_cart_functions[d,d1] = cls.curvilinear_coordinates[d]
                    cls.curvilinear_to_cart_functions[d,d1] = 1
                    cls.curvilinear_to_cart_functions[d,d1] = 1
        return

    def generate_jacobians(cls, coordinate_symbol, order):
        from sympy.tensor.array import MutableDenseNDimArray
        eq = Equation()
        jacobians = "Eq(Jac_i_j, Der(xi_i,x_j))"
        jacobians = eq.expand(jacobians , cls.ndim, coordinate_symbol, substitutions=[], constants=[], Metric=None)
        fdfn = jacobians[:]
        SD = "Eq(sd_i_j_k , Conservative(xi_i,x_j,xi_k))"
        #sdfn = eq.expand(SD , cls.ndim, "x", substitutions=[], constants=[], Metric=None)
        fd_i_j = Matrix(cls.ndim,cls.ndim, lambda i,j:0)
        sd_ijk = MutableDenseNDimArray.zeros(cls.ndim, cls.ndim, cls.ndim)
        sdfn = MutableDenseNDimArray.zeros(cls.ndim, cls.ndim, cls.ndim)
        for i, curv in enumerate(cls.curvilinear_metric):
            for j,stretch in enumerate(cls.stretching_metric):
                if curv and stretch:
                    # IF curvilinear and stretched
                    fd_i_j[i,j] = 1
                    linear_index = i*cls.ndim+j
                    fdfn[linear_index] = Eq(fdfn[linear_index].rhs, fd_i_j[i,j]*fdfn[linear_index].lhs, evaluate=False)
                elif stretch:
                    # if stretched only
                    fd_i_j[j,j] = 1
                    linear_index = i*cls.ndim+j
                    fdfn[linear_index] = Eq(fdfn[linear_index].rhs, fd_i_j[i,j]*fdfn[linear_index].lhs, evaluate=False)
                else:
                    if i==j:
                        linear_index = i*cls.ndim+j
                        fdfn[linear_index] = Eq(fdfn[linear_index].rhs, 1)
                    else:
                        linear_index = i*cls.ndim+j
                        fdfn[linear_index] = Eq(fdfn[linear_index].rhs, 0)
                for k in range(cls.ndim):
                    linear_index = i*cls.ndim+j*cls.ndim + k
                    cur = cls.curvilinear_metric[i]
                    stret = cls.stretching_metric[j]
                    stret1 =cls.curvilinear_metric[k]
                    linear_index = i*cls.ndim*cls.ndim+j*cls.ndim + k
                    if curv and stret and stret1:
                        print i,j,k, "CSS"
                        sd_ijk[i,j,k] = DataObject("SDxi%dx%dxi%d"%(i,j,k))
                        sdfn[i,j,k] = CentralDerivative(DataObject("xi%d"%i),CoordinateObject("x%d"%j) ,DataObject("xi%d"%k))
                        #sdfn[linear_index] = Eq(sdfn[linear_index].lhs, sd_ijk[i,j,k])
                    elif curv and stret:

                        print i,j,k, "CSN"
                        sd_ijk[i,j,k] = DataObject("SDxi%dx%dxi%d"%(i,j,k))
                        sdfn[i,j,k] = CentralDerivative(DataObject("xi%d"%i),CoordinateObject("x%d"%j) ,DataObject("xi%d"%k))
                    elif stret and stret1:
                        print i,j,k, "SS"
                        sd_ijk[k,k,k] = DataObject("SDxi%dx%dxi%d"%(i,j,k))
                        #sdfn[linear_index] = Eq(sdfn[linear_index].lhs, sd_ijk[i,j,k])
                        sdfn[k,k,k] = CentralDerivative(DataObject("xi%d"%i),CoordinateObject("x%d"%j) ,DataObject("xi%d"%k))
                    else:
                        print i,j,k, "NO"
                        sd_ijk[i,j,k] = 0
                        sdfn[i,j,k] = CentralDerivative(DataObject("xi%d"%i),CoordinateObject("x%d"%j) ,DataObject("xi%d"%k))
                        #print sdfn[linear_index], linear_index, i,j,k
                    #sdfn[linear_index] = Eq(sdfn[linear_index].rhs, sd_ijk[i,j,k]*sdfn[linear_index].lhs)
                    #pprint(sdfn[linear_index])

        #sdfn = MutableDenseNDimArray(sdfn,(cls.ndim, cls.ndim, cls.ndim))
        fdfn_subs = cls.convert_to_dictionary(fdfn)
        first_derivative = cls.generate_FD(coordinate_symbol)
        #first_derivative = [f.xreplace(fdfn_subs) for f in first_derivative]
        sd_subs = dict(zip(flatten(sdfn.tolist()), flatten(sd_ijk.tolist())))
        sd = cls.generate_SD(coordinate_symbol)
        #sd = [f.xreplace(fdfn_subs) for f in sd]
        #sd = [f.xreplace(sd_subs) for f in sd]
        #sd = cls.convert_dataobjects(sd)
        #pprint(srepr(sd[0]))
        #pprint(srepr(flatten(sdfn.tolist())[0]))
        #pprint(sd_subs)
        #pprint(sd_ijk)
        import itertools
        # First derivative
        s = [cls.curvilinear_metric, cls.stretching_metric]
        metric_list = list(itertools.product(*s))
        s1 = [[0,1,2], [0,1,2]]
        indices = list(itertools.product(*s1))
        fd_fns = indices[:]
        s1 = [cls.curvilinear_coordinates, cls.cartesian_coordinates]
        fder_fns = (list(itertools.product(*s1)))
        for no,m in enumerate(metric_list):
            #print(der_fns[no], type(der_fns[no]))
            string = ''.join([str(s) for s in flatten(list(fder_fns[no]))])
            fder_fns[no] = CentralDerivative(*fder_fns[no])
            if all(m):
                fd_fns[no] = DataObject("D%s"%string)
            elif m[1] and len(set(indices[no]))==1:
                fd_fns[no] = DataObject("D%s"%string)
            #elif all(not element for element in m):
                #fd_fns[no] = 1
            else:
                fd_fns[no] = 0

        fd_subs = dict(zip(fder_fns, fd_fns))

        # Second derivative
        s = [cls.curvilinear_metric, cls.stretching_metric, cls.curvilinear_metric]
        metric_list = list(itertools.product(*s))
        s1 = [[0,1,2], [0,1,2], [0,1,2]]
        indices = list(itertools.product(*s1))
        sd_fns = indices[:]
        s1 = [cls.curvilinear_coordinates, cls.cartesian_coordinates, cls.curvilinear_coordinates]
        der_fns = (list(itertools.product(*s1)))
        for no,m in enumerate(metric_list):
            #print(der_fns[no], type(der_fns[no]))
            string = ''.join([str(s) for s in flatten(list(der_fns[no]))])
            print m
            der_fns[no] = CentralDerivative(*der_fns[no])
            if all(m):
                sd_fns[no] = DataObject("SD%s"%string)
            elif m[1] and len(set(indices[no]))==1:
                sd_fns[no] = DataObject("SD%s"%string)
            #elif all(not element for element in m):
                #sd_fns[no] = 0
            else:
                sd_fns[no] = 0
        #pprint(der_fns)
        #pprint(sd_fns)
        sd_subs = dict(zip(der_fns, sd_fns))
        #pprint(sd[-1])
        sd = [f.xreplace(fd_subs) for f in sd]
        sd = [f.xreplace(sd_subs) for f in sd]
        #pprint(srepr(der_fns[0]))
        #pprint(first_derivative[0])
        #pprint(srepr(first_derivative[0]))
        #print [(d, d.get_direction) for  d in first_derivative[0].atoms(CentralDerivative)]
        #print [(d, d.direction) for d in cls.curvilinear_coordinates]
        """IT is simple to apply the transformations based on the metric parameters
        what to do is, below everything is a differnt function
        a. generate the differential of the first derivative for an arbitarary function
        b. Genrate the first derivative transformation matrices
        c. apply the first derivative transformation matrices
        d. All the above will be wrapped with the helo of a calling function
        e. A function to transform the given derivative
        f. Similar procedure for the second derivatives
        """
        #exit()
        from .latex import *
        latex = LatexWriter()
        latex.open('./kernels3d.tex')
        metadata = {"title": "Transformations", "author": "Jammy", "institution": ""}
        latex.write_header(metadata)
        latex.write_expression(sd)
        #temp = flatten(sdfn.tolist())
        #temp2 = flatten(sd_ijk.tolist())
        #temp1 = [Eq(x,y, evaluate=False) for x,y in zip(temp2, temp)]
        #latex.write_expression(temp1)
        latex.write_footer()
        latex.close()
        #k = sd_subs.keys()[0]
        #pprint(srepr(k))
        exit()
        return


    def convert_dataobjects(cls, list_of_eq):
        out = []
        #for e in list_of_eq:

        return

    def dot_prodcut(cls, m1, m2):
        out = []
        for no, m in enumerate(m1):
            out += [m*m2[no]]
        out = (type(m1)(out).reshape(*m1.shape))
        return out

    def convert_to_dictionary(cls, list_of_eq):
        dictionary = {}
        for e in list_of_eq:
            dictionary[e.lhs] = e.rhs
        return dictionary

    def generate_FD(cls, coordinate_symbol):
        eq = Equation()
        fd_fn = "Eq(fn_i, MetricDer(f,x_i))"
        fdfn = eq.expand(fd_fn, cls.ndim, coordinate_symbol, substitutions=[], constants=[], Metric=None)
        return fdfn

    def generate_SD(cls, coordinate_symbol):
        eq = Equation()
        sd_fn = "Eq(Der(f, x_i,x_j), MetricDer(MetricDer(f,x_i), x_j))"
        sdfn = eq.expand(sd_fn, cls.ndim, coordinate_symbol, substitutions=[], constants=[], Metric=None)
        return sdfn



    def update_stretching_function(cls, stretch):
        for ind, val in enumerate(stretch):
            if val:
                list1 = cls.cart_to_curvilinear_functions[ind,ind]
                list1[ind] = cls.curvilinear_coordinates[ind]
                list2 = cls.curvilinear_to_cart_functions[ind,ind]
                list2[ind] = cls.cartesian_coordinates[ind]
            else:
                list1 = cls.cart_to_curvilinear_functions[ind,ind]
                list1[ind] = False
                list2 = cls.curvilinear_to_cart_functions[ind,ind]
                list2[ind] = False
        return

    def update_curvilinear_function(cls, curv):
        for ind, val in enumerate(curv):
            if not val:
                for d in range(cls.ndim):
                    list1 = cls.cart_to_curvilinear_functions[ind,d]
                    list2 = cls.cart_to_curvilinear_functions[d,ind]
                    list3 = cls.curvilinear_to_cart_functions[ind,d]
                    list4 = cls.curvilinear_to_cart_functions[d,ind]
                    if d==ind:
                        for dim in range(cls.ndim):
                            list1[dim] = False
                            list2[dim] = False
                            list3[dim] = False
                            list4[dim] = False
                    else:
                        for dim in range(cls.ndim):
                            list1[dim] = False
                            list2[dim] = False
                            list3[dim] = False
                            list4[dim] = False
                    for d1 in range(cls.ndim):
                        list1 = cls.cart_to_curvilinear_functions[d,d1]
                        list1[ind] = False
                        list3 = cls.curvilinear_to_cart_functions[d,d1]
                        list3[ind] = False
        return


    def add_fd_metric_eq(cls, fdmetric):
        cls.FD = fdmetric
        return