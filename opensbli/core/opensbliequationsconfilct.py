from .opensbliobjects import DataSet, ConstantObject
from .opensblifunctions import *
from sympy import flatten, preorder_traversal
from sympy import Equality
class Descritisation(object):
    """
    This should contain the following, TODO descri-->discre
    1. Required_datasets these are the required data sets used in the simulations
    2. Required_constants these are the required constants for the SimulationEquations
    3. Local_evaluations this returns a list of evluations that should be performed for
    the Descritisation of any equations
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


class OpenSBLIEquation(Descritisation):
    pass
class Solution(object):
    def __init__(self):
        self.Kernels = []
        self.datasetranges = {}
        return

class SimulationEquations(OpenSBLIEquation, Solution):
    """
    To proceed
    We will create solution equations, these are nothing but the
    steps required to solve the equations. For example,
    to solve say continuity equation, we already have (rhou0:ndim)
    so the equations would be
    Eq(CD(rhou0,x0), ITS descritised formula)
    similarly, for the momentum equation
    or, the best way would be create an object for evaluating for each and every
    term of the equation CD / WD, TD also includes diagnostic terms
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
    def __new__(cls):
        ret = super(SimulationEquations,cls).__new__(cls)
        ret.equations = []
        return ret

    def add_equations(cls, equation):
        equation = cls._sanitise_equations(equation)
        cls.equations += [equation]
        return

    def descritise(cls, schemes, block):
        """
        Descritisation will be done by
        a. Update the work arrays (based on location) for the CD, WD, TD --> irrespective of Type
        b. Set the evaluation range of each of the workbases to be that of the grid, irrespective of CD, WD or TD
        c. Descritise each and every CD, WD, TD
        """
        (Solution,cls).__init__(cls)
        #all_discretisations = cls.required_functions_local
        #for d in all_discretisations:
            #if not block.store_derivatives:
                #d.store = False
            #elif d in block.derivatives_to_store:
                #d.store = False
            #else:
                #d.update_work(block)
        # Required DataSetBases


        # Descritise the derivatives
        for sc in schemes:
            schemes[sc].discretise(cls, block)
        ## Now update the range of evaluations of the DataSetBases
        #dbases = []
        #for d in all_discretisations:
            #dbases += d.required_datasetbases
        ## set the range of evaluation o
        ##pprint(dbases)
        #pprint(set(dbases))
        #for sc in schemes:
            #pprint(schemes[sc].scheme_required_databases)
            #pprint(schemes[sc].equations)
        return
"""
1. Convert CD(u0,x0,x1)--> CD(CD(u0,x0),x1) :: This is applied for all functions that are not homogeneous
2. Do a tree traversal of expressions
    a. if isinstance(expr,DataSet):
        pass
    b. elif do post traverse of tree for any CD
    c. if CD found apply the CD, and create a evaluation expression
    d. get the expression without CD and apply the outer CD
3. The simple way to do all these is to perform a post order tree traversal and findout CDs and create equaitons for
   inner CD and finally outer CD
4. While doing this tree traversal, update dataset range of evaluation.
This completes CD evaluation.

For TD evaluation, follow  a similiar procedure
For WD evaluation,
1. If it is vector type reconstruction, pass in all equations and do the procedure.
2, If scaler type, do a procedure similar lo CD

Finally, update the original equaitons with their respective work DataSets

Using the same classes diagnostics metric equaitons or any other can be applied which then can be passed to algorith to generate the final kernels which inturn can be used to writeout the final code.


Once you do this, you are on top of the world.
"""

class DiagnosticsEquation(OpenSBLIEquation):

    def __new__(cls):

        return

class MetricsEquation(OpenSBLIEquation):

    def __new__(cls):

        return


