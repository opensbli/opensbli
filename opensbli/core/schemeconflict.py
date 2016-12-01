
from sympy.calculus import finite_diff_weights
from sympy import *
from sympy.printing import *
from .opensbliobjects import ConstantObject
from .opensblifunctions import CentralDerivative
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

class Central(Scheme):

    """ Spatial discretisation scheme using central differences. """

    def __init__(self, order):
        """ Set up the scheme.

        :arg int order: The order of accuracy of the scheme.
        """
        Scheme.__init__(self, "CentralDerivative", order)
        # Points for the spatial scheme
        self.points = list(i for i in range(-order/2, order/2+1))
        # Set max_order to 2 currently
        max_order = 2
        #self._generate_derivative_weights(max_order)
        self.required_database = []
        return

    def _generate_weights(self, direction, order):
        """ Descritises only the homogeneous derivatives of any order or
        first derivatives"""
        #print(self.points, direction, type(direction))
        self.diffpoints = [i*ConstantObject('Delta_%s'%(direction)) for i in self.points]
        weights = finite_diff_weights(order, self.diffpoints, 0)
        return weights[order][-1]

    def add_required_database(self, dbases):
        self.required_database += flatten(list(dbases))
        return
    @property
    def scheme_required_databases(self):
        return set(self.required_database)
    def update_works(self, to_descritse, block):
        #for d in to_descritise:
            #if not block.store_derivatives:
                    #d.store = False
                #elif d in block.derivatives_to_store:
                    #d.store = False
                #else:
                    #d.update_work(block)
        # How to do traditional CFD
        return
    def discretise(self, type_of_eq, block):
        """
        Do here the following
        a. Find all the Function of type(self.name)
        b. Add all the DataSetBases required to the type_of_eq
        c. Create Equations for the evaluation ro create Kernels of each function depending on the grid/block
            control parameters
        d. Set the range of evaluation of the DataSetBases
        e. Update the Descritised equations in type_of_eq by substituting the equations with respective
            work array or discretised formula
        """
        fns = type_of_eq.required_functions_local
        required_central_ders = []
        for fn in fns:
            if isinstance(fn, CentralDerivative):
                required_central_ders += [fn]
        to_discretise = set(required_central_ders)
        pprint(to_discretise)
        #if not
        #for d in to_discretise:
            #if isinstance(d, CentralDerivative):
                #expr = self.traverse(d)
                ## update the expression for the discretisation
                #d.discretised_expr = expr._discretise_derivative(self)
                #self.add_required_database(d.required_datasetbases)
                ##self.add_required_database(expr.required_datasetbases)
                #self.equations += [Eq(d.work, d.discretised_expr)]
        return

    def traverse(self, CD):
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
            for cd in inner_cds[:-1]:
                if cd.work:
                    self.add_required_database(cd.required_datasetbases)
                    expr = expr.subs(cd, cd.work)
                    #self.add_required_database(cd.required_datasetbases)
                else:
                    # THIS raises an error when the CD(u0,x0) is not there in all derivatives ,
                    # while evaluating CD(CD(u0,x0),x1)
                    raise ValueError("NOT IMPLEMENTED THIS")

        return expr


    def set_halos(self, block):
        halos = (-self.order/2,self.order/2)
        self.halos = block.set_block_halos(halos)
        # Copy the block indices, mapping and deltax_i as they are reused frequently
        self.derivative_direction = block.indices
        self.index_mapping = block.mapped_indices
        self.deltas = block.deltas
        return

    def setup(self, block):
        self.set_halos(block)
        equations = flatten(self.scheme_equations)
        spatial, temporal = get_derivatives(equations)
        self.local_derivatives = spatial
        self.indexed_required =  list(set(flatten([list(eq.atoms(Indexed)) for eq in equations])))
        return

