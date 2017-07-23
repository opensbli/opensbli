from opensbli.core.parsing import Equation
from opensbli.core.latex import *
from opensbli.core.opensbliequations import NonSimulationEquations,Discretisation, Solution, OpenSBLIEquation
from opensbli.core.opensblifunctions import CentralDerivative
from opensbli.initialisation.common import *
from sympy.tensor.array import MutableDenseNDimArray
from opensbli.core.opensbliobjects import *
from opensbli.core.kernel import Kernel

def convert_dataset_base_expr_to_datasets(expression, index):
        for a in expression.atoms(DataSet):
            b = a.base
            expression = expression.xreplace({a: b[index]})
        return expression

class MetricsEquation(NonSimulationEquations,Discretisation, Solution):
    def __new__(cls, **kwargs):
        ret = super(MetricsEquation,cls).__new__(cls, **kwargs)
        ret.equations = []
        ret.kwargs = {'strong_differentiability':True}
        ret.algorithm_place = [BeforeSimulationStarts()]
        ret.order = 1
        return ret
    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h
    def _hashable_content(self):
        return "MetricsEquation"

    def genreate_transformations(cls, ndim, coordinate_symbol, parameters, max_order):
        cls.ndim = ndim
        if len(flatten(parameters)) != ndim*2:
            raise ValueError("The parameters for stretching provided should match the number of dimensions")
        cls.transformation_eq = [0 for i in range(max_order)]
        cls.stretching_metric = [param[0] for param in parameters] # Get whether true or false for streching
        cls.curvilinear_metric = [param[1] for param in parameters] # Get T/F for curvilinear
        cart = CoordinateObject('%s_i'%(coordinate_symbol))
        curv = CoordinateObject('%s_i'%('xi'))
        cartesian_coordinates = [cart.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        curvilinear_coordinates = [curv.apply_index(cart.indices[0], dim) for dim in range(ndim)]
        cls.curvilinear_coordinates = curvilinear_coordinates
        cls.cartesian_coordinates = cartesian_coordinates
        cls.latex_debug_start()
        fd_subs, fd_fns = cls.transform_first_derivative(coordinate_symbol)
        cls.transform_second_derivative(coordinate_symbol, fd_subs, fd_fns)
        cls.latex_debug_end()
        return

    def latex_debug_start(self):
        self.latex_file = LatexWriter()
        latex = self.latex_file
        latex.open('./metric_transformations.tex')
        metadata = {"title": "Transformations of the equations in OpenSBLI framework", "author": "Satya P Jammy", "institution": "University of Southampton"}
        latex.write_header(metadata)
        latex.write_string("The Cartesian coordinates system is (%s)"%(','.join([str(s) for s in self.cartesian_coordinates])))
        latex.write_string("The Curvilinear coordinate system is vector is (%s)"%(','.join([str(s) for s in self.curvilinear_coordinates])))
        if any(self.curvilinear_metric):
            latex.write_string("The transformation parameters provided are %s"%(','.join(["curvilinear in "+ str(self.cartesian_coordinates[no]) for no,s in enumerate(self.curvilinear_metric) if s])))
        latex.write_string("The transformation parameters provided are %s"%(','.join(["stretched in "+ str(self.cartesian_coordinates[no]) for no,s in enumerate(self.stretching_metric) if s])))
        return

    def latex_debug_end(self):
        latex = self.latex_file
        latex.write_footer()
        latex.close()
        return

    def transform_first_derivative(cls, coordinate_symbol):
        A = Matrix(cls.ndim, cls.ndim,lambda i,j:cls.cartesian_coordinates[i])
        M1 = Matrix(cls.ndim, cls.ndim,lambda i,j:cls.curvilinear_coordinates[j])
        M2 = Matrix(cls.ndim, cls.ndim,lambda i,j:0)
        # Full 3D curvilinear expansion of the first derivatives
        fds = cls.full3D_FD_transformation(coordinate_symbol)
        cls.detJ = DataObject("detJ")
        fd_subs = {}
        fd_symbols = {}
        fd_metric_eqs = []
        eq = Equation()
        fd_jacobians = Matrix(cls.ndim, cls.ndim,lambda i,j:DataObject('D%d%d'%(i,j)))
        Cartesian_curvilinear_derivatives  = zeros(cls.ndim, cls.ndim)

        # Apply the curvilinear metrics
        for i, curv in enumerate(cls.curvilinear_metric):
            if curv:
                for j, curv2 in enumerate(cls.curvilinear_metric):
                    if curv and curv2:
                        M2[i,j] = cls.cartesian_coordinates[j]

        # update the stretching metrics
        for i, curv in enumerate(cls.curvilinear_metric):
            if not curv:
                for j, stretch in enumerate(cls.stretching_metric):
                    if i ==j:
                        if stretch:
                            M2[i,j] = cls.cartesian_coordinates[j]
                        else:
                            M2[i,j] = cls.curvilinear_coordinates[j]
                    else:
                        M2[i,j] = 0
        for i in range(cls.ndim):
            for j in range(cls.ndim):
                args_orig = [cls.curvilinear_coordinates[i],cls.cartesian_coordinates[j]]
                args_eval = [cls.curvilinear_coordinates[i],M2[i,j]]
                args_metricder = [DataObject("%s"%cls.cartesian_coordinates[i]), cls.curvilinear_coordinates[j]]
                v = CentralDerivative(*args_eval).doit()
                if isinstance(v,CentralDerivative):
                    fd_subs[CentralDerivative(*args_orig)] =  fd_jacobians[i,j]
                    Cartesian_curvilinear_derivatives[i,j] = CentralDerivative(*args_metricder)
                else:
                    fd_subs[CentralDerivative(*args_orig)] = v
                    fd_jacobians[i,j] = v
                    Cartesian_curvilinear_derivatives[i,j] = v
        fd_transformed = zeros(cls.ndim,1)
        for f in fds:
            d = (f.lhs.get_direction)[0]
            f = f.subs(fd_subs)
            fd_transformed[d] = f.rhs

        cls.FD_transformation_equations = fd_transformed
        cls.transformation_eq[0] = fd_transformed
        cls.FD_metrics = fd_jacobians
        cls.generate_fd_metrics_equations(Cartesian_curvilinear_derivatives)

        cls.classical_strong_differentiabilty_transformation = []
        for d in fd_transformed:
            # Donot convert the convective derivatives
            #cls.classical_strong_differentiabilty_transformation += [cls.convert_der(d)]
            cls.classical_strong_differentiabilty_transformation += [d]
        # Write latex file for easy debugging
        latex = cls.latex_file
        for i in range(cls.ndim):
            cd = CentralDerivative(cls.general_function, cls.cartesian_coordinates[i])
            latex.write_expression(Eq(cd, cls.classical_strong_differentiabilty_transformation[i]))
        #latex.write_expression(cls.FD_metrics)
        #latex.write_expression(Cartesian_curvilinear_derivatives.adjugate())
        #latex.write_expression(M2)
        return fd_subs, M2

    def convert_der(self, expr):
        """ Converts any function v*D(f,x) to D(v*f,x) - f*D(v,x)
        Below substitution, returns equation F-31 of http://cfl3d.larc.nasa.gov/Cfl3dv6/V5Manual/GenCoor.pdf
        >>> expr,matches = expr.replace(v*Derivative(f,x), Derivative(v*f,x) -f*Derivative(v,x),map=True,simultaneous=False, exact=True)
        As we know that -f*D(v,x) gets cancelled out where v is (xi_ij/J), We return only
        D(v*f,x)
        """
        v1 = Wild('x')
        f = WildFunction('F', nargs = 2)
        # Get all the functions with arguments 2, i.e first derivatives
        l = list(expr.find(f))
        # find all the atoms in the expression i.e a DataObject and S.One as some of the terms have no metric
        l1 = [a for a in list(expr.find(v1)) if isinstance(a,DataObject)] + [S.One]
        terms = [a*b for a in l for b in l1]
        # This is D(v*f,x) - f*D(v,x)
        #term_subs = [(type(self)(a.args[0]*b, a.args[1:]) -a.args[0]*type(self)(b, a.args[1:]))  for a in l for b in l1]
        # As f*D(v,x) get cancelled out in the summation return only v*D(f,x)
        #term_subs = [CentralDerivative(a.args[0]*b/self.detJ , a.args[1:])  for a in l for b in l1]
        #expr = self.detJ*expr.subs(dict(zip(terms, term_subs))) # substitute
        term_subs = [CentralDerivative(a.args[0]*b*self.detJ , a.args[1:])  for a in l for b in l1]
        pprint(expr)
        expr = expr.subs(dict(zip(terms, term_subs)))/self.detJ # substitute
        pprint(expr)
        # WARNING check this expression
        return expr


    def generate_fd_metrics_equations(cls, Cartesian_curvilinear_derivatives):
        adjointJ = Cartesian_curvilinear_derivatives.adjugate(); detJ = Cartesian_curvilinear_derivatives.det()
        evaluation = adjointJ/detJ
        eqns = []
        # TODO add if condition on the arguments if this is required
        eqns = [Eq(cls.detJ, detJ)]
        eqns += [Eq(x,y) for (x,y) in zip(cls.FD_metrics, evaluation)]
        eqns = [e for e in eqns if isinstance(e, Equality)]
        eqns = [OpenSBLIEquation(eq.lhs, eq.rhs) for eq in eqns]
        cls.fdequations = eqns
        return eqns

    def transform_second_derivative(cls, coordinate_symbol, fd_subs, M2):
        # Second derivative substitutions
        SD = MutableDenseNDimArray([0 for i in range(cls.ndim) for j in range(cls.ndim) for k in range(cls.ndim)],
                   (cls.ndim, cls.ndim, cls.ndim))
        sd_subs = {}
        SD_evaluations =  MutableDenseNDimArray([0 for i in range(cls.ndim) for j in range(cls.ndim) for k in range(cls.ndim)], (cls.ndim, cls.ndim, cls.ndim))
        SD_jacobians = MutableDenseNDimArray([0 for i in range(cls.ndim) for j in range(cls.ndim) for k in range(cls.ndim)],
                   (cls.ndim, cls.ndim, cls.ndim))
        latex = cls.latex_file
        for i in range(M2.shape[0]):
            for j in range(M2.shape[0]):
                #if M2[i,j] == 0:
                    #pass
                if M2[i,j] in cls.curvilinear_coordinates:
                    pass
                else:
                    # Apply the metrics
                    for k, curv in enumerate(cls.curvilinear_metric):
                        if curv:
                            SD[i,j,k] = cls.cartesian_coordinates[k]
                        elif cls.stretching_metric[k]:
                            if i==j==k:
                                SD[i,j,k] = cls.cartesian_coordinates[k]
                for k in range(cls.ndim):
                    args_orig = [cls.curvilinear_coordinates[i],cls.cartesian_coordinates[j], cls.curvilinear_coordinates[k]]
                    args_eval = [cls.curvilinear_coordinates[i],cls.cartesian_coordinates[j], SD[i,j,k]]
                    v = CentralDerivative(*args_eval).doit()
                    if isinstance(v,CentralDerivative):
                        sd_subs[CentralDerivative(*args_orig)] = DataObject("SD%d%d%d"%(i,j,k))
                        SD_jacobians[i,j,k] = DataObject("SD%d%d%d"%(i,j,k))
                        derivative_args = [cls.FD_metrics[i,j], cls.curvilinear_coordinates[k]]
                        SD_evaluations[i,j,k] = CentralDerivative(*derivative_args)
                    else:
                        sd_subs[CentralDerivative(*args_orig)] = v
                        SD_jacobians[i,j,k] = v
                        SD_evaluations[i,j,k] = v

        # Full 3D curvilinear expansion of the second derivatives
        sds = cls.full3D_SD_transformation(coordinate_symbol)
        sd_transformed = zeros(cls.ndim,cls.ndim)
        for s in sds:
            d = (s.lhs.get_direction)
            s = s.subs(fd_subs).subs(sd_subs)
            sd_transformed[tuple(d)] = s.rhs
        cls.SD_transformation_equations = sd_transformed
        cls.classical_strong_differentiabilty_transformation_sd = zeros(cls.ndim, cls.ndim)
        for no,d in enumerate(sd_transformed):
            # if it is a second derivative we will just multiply by the detJ as of now
            cls.classical_strong_differentiabilty_transformation_sd[no] = d

        cls.SD_metrics = SD_jacobians
        cls.SD_evaluations = SD_evaluations
        cls.generate_sd_metrics_equations()
        #fn_ders = CentralDerivative(DataObject('f'), cls.cartesian_coordinates[i], cls.cartesian_coordinates[j]) for i,j in zip(()
        for i in range(cls.ndim):
            for j in range(cls.ndim):
                fn = CentralDerivative(cls.general_function, cls.cartesian_coordinates[i], cls.cartesian_coordinates[j])
                latex.write_expression(Eq(fn,cls.classical_strong_differentiabilty_transformation_sd[i,j]))

        #eqns = [Eq(x,y,evaluate=False) for x,y in zip(cls.SD_metrics, cls.SD_evaluations)]
        #for eq in eqns:
            #latex.write_expression(eq)
        #pprint(srepr(cls.SD_transformation_equations[0]))
        return
    def generate_sd_metrics_equations(cls):
        eqns = [Eq(x,y) for x,y in zip(cls.SD_metrics, cls.SD_evaluations)]
        eqns = [e for e in eqns if isinstance(e, Equality)]
        eqns = [OpenSBLIEquation(eq.lhs, eq.rhs) for eq in eqns]
        cls.sdequations = eqns
        return

    def full3D_FD_transformation(cls,coordinate_symbol):
        eq = Equation()
        fd_fn = "Eq(Der(f, %s_i), MetricDer(f,%s_i))"%(coordinate_symbol,coordinate_symbol)
        cls.general_function = CoordinateObject('f')
        fdfn = eq.expand(fd_fn, cls.ndim, coordinate_symbol, substitutions=[], constants=[], Metric=None)
        return fdfn

    def full3D_SD_transformation(cls, coordinate_symbol):
        eq = Equation()
        sd_fn = "Eq(Der(f, %s_i,%s_j), MetricDer(MetricDer(f,%s_i), %s_j))"%(coordinate_symbol,coordinate_symbol,coordinate_symbol,coordinate_symbol)
        sdfn = eq.expand(sd_fn, cls.ndim, coordinate_symbol, substitutions=[], constants=[], Metric=None)
        return sdfn

    def apply_transformation(cls, equation):
        if cls.kwargs['strong_differentiability']:
            fns = equation.atoms(Function)
            #pprint(fns)
            for fn in equation.atoms(Function):
                transformed = fn.classical_strong_differentiabilty_transformation(cls)
                equation = equation.subs({fn:transformed})
        else:
            raise NotImplementedError("Not implemented other transformation of Metric ders")
        return equation


    def create_residual_arrays(cls):
        for eq in flatten(cls.equations):
            eq.residual = eq.lhs
        return

    def spatial_discretisation(cls, schemes, block):
        (Solution,cls).__init__(cls)
        if any(cls.curvilinear_metric):
            raise NotImplementedError("Handling curvilinear coordinates for the evaluation of metrics is not applied")

        #cls.descritsed_equations = copy.copy(cls.equations)
        spatialschemes = []
        for sc in schemes:
            if schemes[sc].schemetype == "Spatial":
                spatialschemes += [sc]
        cls.requires = {}
        block.store_work_index # store the work array index
        for no,sc in enumerate(spatialschemes):
            # TODO discretise the First derivatives and then apply BC's then discretise Second derivatives
            cls.equations  = block.dataobjects_to_datasets_on_block(cls.fdequations) # First derivatives
            cls.create_residual_arrays()
            schemes[sc].required_constituent_relations = {}
            a = schemes[sc].discretise(cls, block)
            #pprint(cls.sdequations)
            cls.apply_periodic_bc(block)
            cls.equations  = block.dataobjects_to_datasets_on_block(cls.sdequations) # Second derivatives
            cls.create_residual_arrays()
            schemes[sc].required_constituent_relations = {}
            b = schemes[sc].discretise(cls, block)
            schemes[sc].required_constituent_relations = {}
            #if b:
                #for key, val in a.iteritems():
                    #cls.requires[key] = val
        for k in cls.Kernels:
            print k.computation_name
        #exit()
        #pprint(cls.requires)
        #pprint([k.equations for k in cls.Kernels])
        #exit()
        cls.process_kernels(block)
        block.reset_work_to_stored # reset the work array index on the block
        #for k in cls.Kernels:
            #print k.computation_name, k.halo_ranges, k.ranges
        #exit()
        return

    def apply_periodic_bc(self, block):
        arrays = self.FD_metrics[:] + [self.detJ]
        arrays = block.dataobjects_to_datasets_on_block(arrays)
        from opensbli.core.bcs import PeriodicBoundaryConditionBlock as pbc
        modify = {}
        for no, val in enumerate(block.boundary_types):
            left = val[0]; right = val[1]
            if isinstance(left, pbc):
                if no in modify:
                    modify[no][0] = left
                else:
                    modify[no] = [left, None]
            if isinstance(right, pbc):
                if no in modify:
                    modify[no][1] = right
                else:
                    modify[no] = [None, right]

        if modify:
            for dire in modify:
                boundary_mods = [k for k in modify[dire] if k] 
                for b in boundary_mods:
                    self.Kernels += [b.apply(arrays, b.direction, b.side, block)]

        return
    def apply_boundary_conditions(cls, block):
        """No global boundary conditions for the metric terms"""
        ndim = block.ndim
        directions = [i for i in range(ndim)]
        metrics = Matrix(ndim, ndim, cls.FD_metrics[:])
        detJ = cls.detJ
        kernels = []
        # arrays = metrics[:] + [detJ]
        # arrays = block.dataobjects_to_datasets_on_block(arrays)
        # pprint(arrays)
        # kernels = block.apply_boundary_conditions(arrays)
        #for direction in directions:
            #for side in [0,1]:
                #kernels += [create_metric_bc_kernel(metrics, detJ, direction, side, block)]
        #cls.Kernels += kernels
        #exit()

        return

    def process_kernels(cls, block):
        for key,kernel in cls.constituent_relations_kernels.iteritems():
            if isinstance(kernel, Kernel):
                kernel.update_block_datasets(block)
        for kernel in cls.Kernels:
            if isinstance(kernel, Kernel):
                kernel.update_block_datasets(block)
        return


def create_metric_bc_kernel(metrics, detJ, direction, side, block):
    # if block.ndim == 1:
    #     transformed_metrics = -1*metrics
    #     transformed_det = -1*detJ
    # else:
    #     transformation_matrix = eye(block.ndim)
    #     transformation_matrix[direction, direction] = -1
    #     transformed_metrics = metrics*transformation_matrix.adjugate()
    #     transformed_det = transformation_matrix.det()*detJ

    kernel = Kernel(block, computation_name="Metric boundary dir%d side%d" % (direction, side))
    kernel.set_boundary_plane_range(block, direction, side)
    kernel.ranges[direction][side] = block.ranges[direction][side]
    halos = kernel.get_plane_halos(block)
    # Add the halos to the kernel in directions not equal to boundary direction
    for i in [x for x in range(block.ndim) if x != direction]:
        kernel.halo_ranges[i][0] = block.boundary_halos[i][0]
        kernel.halo_ranges[i][1] = block.boundary_halos[i][1]
    # # Add halos across the boundary side only
    # kernel.halo_ranges[direction][side] = block.boundary_halos[direction][side]

    loc = [0 for i in range(block.ndim)]
    new_loc = loc[:]
    if side == 0:
        from_side_factor = -1
        to_side_factor = 1
    elif side == 1:
        from_side_factor = 1
        to_side_factor = -1

    lhs = [x for x in metrics if x not in [-1, 0, 1]] + [detJ]
    rhs = lhs[:]
    # rhs = [x for x in transformed_metrics if x not in [-1, 0, 1]] + [transformed_det]
    lhs = block.dataobjects_to_datasets_on_block(lhs)
    rhs = block.dataobjects_to_datasets_on_block(rhs)

    new_loc[direction] += to_side_factor

    transfer_indices = [tuple([from_side_factor*t, to_side_factor*t]) for t in range(1, abs(halos[direction][side]) + 1)]

    final_equations = []
    for index in transfer_indices:
        array_equations = []
        loc_lhs, loc_rhs = loc[:], loc[:]
        loc_lhs[direction] += index[0]
        loc_rhs[direction] += index[1]
        for left, right in zip(lhs, rhs):
            left = convert_dataset_base_expr_to_datasets(left, loc_lhs)
            right = convert_dataset_base_expr_to_datasets(right, loc_rhs)
            array_equations += [Eq(left, right, evaluate=False)]
        final_equations += array_equations
    kernel.add_equation(final_equations)
    pprint(kernel.equations)
    kernel.update_block_datasets(block)
    return kernel

"""
TODO remove the process kernels in all the equation classes
this will be done in algorithm
b. As there is no curvilinear system
"""
