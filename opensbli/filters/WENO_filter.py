""" David J. Lusher 08/2020: WENO non-linear filter for shock-capturing."""
from opensbli import *
from sympy import symbols, exp, pprint, Piecewise, binomial, Min, sqrt, Equality, tanh
from opensbli.core.opensbliobjects import DataObject, ConstantObject, GroupedPiecewise
from opensbli.equation_types.opensbliequations import OpenSBLIEquation, NonSimulationEquations, ConstituentRelations
from opensbli.postprocess.post_process_eq import *
from opensbli.core.kernel import ConstantsToDeclare as CTD
from opensbli.code_generation.algorithm.common import *
from opensbli.utilities.user_defined_kernels import UserDefinedEquations
from opensbli.schemes.spatial.weno import *
from opensbli.core.boundary_conditions.bc_core import WallBC
from sympy.functions.elementary.piecewise import ExprCondPair


class WENOFilter(NonSimulationEquations):
    """ Class to apply a WENO-based non-linear filter after a full time-step of a non-dissipative high order base scheme. The dissipative
    portion of a WENO procedure is used in characteristic space, by substracting a central difference flux approximation of order n+1. The shock location sensor
    uses the absolute difference of the non-linear to ideal WENO weights. The amount of dissipation is controlled by Mach number or dilatation/vorticity sensors. The governing
    equations in the user script should be central derivatives in a skew-symmetric formulation to improve numerical stability."""
    def __init__(self, block, order, metrics=None, dissipation_sensor='Ducros', Mach_correction=False):
        self.reconstruction_kernels = []
        self.residual_kernels = []
        block.shock_filter = True
        # Choice of how to evaluate the amount of dissipation to be added (varies spatially in the domain)
        self.dissipation_sensor = dissipation_sensor
        self.Mach_correction = Mach_correction
        self.block = block
        self.ndim = block.ndim
        self.order = order
        self.equation_classes = []
        self.nhalos = 5
        # Store the value of the dissipation control sensor for debugging
        self.store_kappa = True
        # Scheme used to form the non-linear filter
        self.scheme_type = "**{\'scheme\':\'Weno\'}"
        # Check if the problem needs a metric transformation of the equations
        self.process_metrics(metrics)
        self.constants = ["Re", "Pr","gama", "Minf", "SuthT", "RefT"]
        # Ensure gama has been added to the constants to define
        self.gama = ConstantObject('gama')
        CTD.add_constant(self.gama)
        # Einstein class to expand equations
        self.EE = EinsteinEquation()
        # Main class to generate the filter
        self.main(order, block)
        return

    def detect_wall_boundaries(self):
        """ The shock-filter is turned off in the near-wall region. This function detects which directions, if any, have
        wall boundary conditions."""
        self.wall_boundaries = [[False, False] for _ in range(self.ndim)]
        try:
            for direction in range(self.ndim):
                for side in [0,1]:
                    if isinstance(self.block.boundary_types[direction][side], WallBC):
                        self.wall_boundaries[direction][side] = True
        except:
            raise ValueError("Please set boundary conditions on the block before calling the shock filter.")
        return

    def process_metrics(self, metrics):
        # Uniform mesh
        if metrics is None:
            self.metric_class = None
            self.stretched, self.curvilinear = False, False
        # Stretched or curvlinear mesh
        else:
            self.metric_class = metrics
            # Check whether the mesh is only stretched or full curvilinear
            if sum(self.metric_class.stretching_metric) > 0:
                self.stretched = True
            if sum(self.metric_class.curvilinear_metric) > 0:
                self.curvilinear = True
            else:
                self.curvilinear = False
        return

    def create_weno_equations(self):
        # Define the compresible Navier-Stokes equations in Einstein notation, depending on the metric input
        scheme_type = self.scheme_type
        # Uniform mesh, no stretching or curvilinear terms
        if self.metric_class is None:
            coordinate_symbol = "x"
            mass = "Eq(Der(rho,t), - Conservative(rho*u_j,x_j,%s))" % scheme_type
            momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s))" % scheme_type
            energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s))" % scheme_type
            output_equations = flatten([self.EE.expand(eq, self.ndim, coordinate_symbol, [], self.constants) for eq in flatten([mass, momentum, energy])])
        else:
            # Full curvilinear
            if self.curvilinear:
                coordinate_symbol = "xi"
                a = "Conservative(detJ * rho*U_j,xi_j,%s)" % scheme_type
                mass = "Eq(Der(rho,t), - %s/detJ)" % (a)
                a = "Conservative(detJ * (rhou_i*U_j + p*D_j_i), xi_j , %s)" % scheme_type
                momentum = "Eq(Der(rhou_i,t) , -  %s/detJ)" % (a)
                a = "Conservative(detJ * (p+rhoE)*U_j,xi_j, %s)" % scheme_type
                energy = "Eq(Der(rhoE,t), - %s/detJ)" % (a)
                output_equations = flatten([self.EE.expand(eq, self.ndim, coordinate_symbol, [], self.constants) for eq in flatten([mass, momentum, energy])])
            # Only stretching is applied
            else:
                coordinate_symbol = "x"
                mass = "Eq(Der(rho,t), - Conservative(rho*u_j,x_j,%s))" % scheme_type
                momentum = "Eq(Der(rhou_i,t) , -Conservative(rhou_i*u_j + KD(_i,_j)*p,x_j , %s))" % scheme_type
                energy = "Eq(Der(rhoE,t), - Conservative((p+rhoE)*u_j,x_j, %s))" % scheme_type
                governing_eq = flatten([self.EE.expand(eq, self.ndim, coordinate_symbol, [], self.constants) for eq in flatten([mass, momentum, energy])])
                output_equations = flatten([self.metric_class.apply_transformation(eqn) for eqn in (governing_eq)])                          
        return output_equations

    def create_kernel(self, name, equations, halo_type, block):
        filter_class = UserDefinedEquations()
        filter_class.algorithm_place = InTheSimulation(frequency=False)
        filter_class.computation_name = name
        filter_class.order = self.component_counter
        # Add the halo type to extend the range of evaluation
        filter_class.halos = halo_type
        filter_class.add_equations(equations)
        return filter_class

    def add_kernel(self, kernel):
        """ Adds the finished kernels to the storage."""
        if isinstance(kernel, list):
            for ker in kernel:
                self.equation_classes.append(ker)
        else:
            self.equation_classes.append(kernel)
        return

    def group_by_direction(self, eqs):
        """ Groups the input equations by the direction (x0, x1, ...) they depend upon.

        :arg list eqs: List of equations to group by direction.
        :returns: dict: grouped: Dictionary of {direction: equations} key, value pairs for equations grouped by direction."""
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(WenoDerivative))
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        return grouped

    def constituent_relations(self, block):
        """ Evalutes the constiteunt relations on the state at the end of a full step
        of the Runge-Kutta explicit time-stepper. Only the invscid terms are evaluted here (no viscosity relation)"""
        CR_eqns = []
        # Conservative Q array entries from the current state
        rho, rhoE = self.solution_vector[0], self.solution_vector[-1]
        momentum_components = [self.solution_vector[i+1] for i in range(block.ndim)]
        inv_rho = GridVariable('inv_rho')
        CR_eqns += [OpenSBLIEq(inv_rho, 1.0/rho)]
        # Primitive components and speed of sound
        velocity_components = [block.location_dataset('u%d' % i ) for i in range(block.ndim)]
        CR_eqns += [OpenSBLIEq(x, y*inv_rho) for (x, y) in zip(velocity_components, momentum_components)]
        # Pressure equation
        p, a = block.location_dataset('p'), block.location_dataset('a')
        # rhoE = p/(gama-1) + 0.5*(rhou**2)/rho
        CR_eqns += [OpenSBLIEq(p, (self.gama-1)*(rhoE - 0.5*sum([dset**2 for dset in momentum_components])*inv_rho))]
        CR_eqns += [OpenSBLIEq(a, sqrt(self.gama*p*inv_rho))]
        # Projected velocities if full curvilinear coordinates are being used
        if self.curvilinear:
            metric_vel = "Eq(U_i, D_i_j*u_j)"
            CR_eqns += flatten([self.EE.expand(metric_vel, self.ndim, "x", [], self.constants)])

        # Optiional Low Mach number correction
        if self.Mach_correction:
            M_var, M_eqns = self.evaluate_Yee_Mach_sensor(block)
            CR_eqns += M_eqns

        CR_halos = []
        for _ in range(block.ndim):
            CR_halos.append([self.halo_type, self.halo_type])
        CR_kernel = self.create_kernel('Constituent Relations evaluation', CR_eqns, CR_halos, block)
        self.component_counter += 1
        # Add the kernel
        self.add_kernel(CR_kernel)
        return

    def evaluate_Yee_Mach_sensor(self, block):
        """ Sensor controlling the amount of dissipation to apply. Turns the filter off in low-Mach regions.
        (High Order Filter Methods for Wide Range of Compressible Flow Speeds, Yee, 2010). Primitive variables are computed locally
        to reduce the number of input/output arrays to the kernel."""
        Mach_equations = []
        # Evaluate the local Mach number
        M, inv_a, p, inv_rho = symbols('M inv_a p inv_rho', **{'cls' : GridVariable})
        # Conservative quantities
        rho, rhoE = self.solution_vector[0], self.solution_vector[-1]
        momentum_components = [self.solution_vector[i+1] for i in range(block.ndim)]
        velocity_components = [GridVariable('u%d' % i ) for i in range(block.ndim)]
        Mach_equations += [OpenSBLIEq(inv_rho, 1.0/rho)]
        Mach_equations += [OpenSBLIEq(x, y*inv_rho) for (x, y) in zip(velocity_components, momentum_components)]
        Mach_equations += [OpenSBLIEq(p, (self.gama-1)*(rhoE - 0.5*sum([dset**2 for dset in momentum_components])*inv_rho))]
        Mach_equations += [OpenSBLIEq(inv_a, 1.0/sqrt(self.gama*p*inv_rho))]
        Mach_equations += [OpenSBLIEq(M, sqrt(sum(dset**2 for dset in velocity_components))*inv_a)]
        # Evaluation of the kappa parameter to control the amount of dissipaton
        Mach_correct = block.location_dataset('Mach_sensor')
        Mach_equations += [OpenSBLIEq(Mach_correct, Min(0.5*M**2 * sqrt(4+(1-M**2)**2) / (1+M**2), 1))]
        return Mach_correct, Mach_equations

    def evaluate_Ducros_sensor(self, block):
        # Add a shock sensor for the WENO filter
        SS = ShockSensor()
        output_eqns, kappa = SS.ducros_equations(block, "x", metrics=self.metric_class, name='kappa')
        # Low Mach number correction
        if self.Mach_correction:
            sensor_evaluation = output_eqns[-1]
            del output_eqns[-1]
            output_eqns += [OpenSBLIEq(sensor_evaluation.lhs, block.location_dataset('Mach_sensor')*sensor_evaluation.rhs)]

        # Halo points for the sensor kernel
        sensor_halos = []
        for _ in range(block.ndim):
            sensor_halos.append([self.halo_type, self.halo_type])
        sensor_kernel = self.create_kernel('Shock sensor', flatten(output_eqns), sensor_halos, block)
        # Add the kernel
        self.add_kernel(sensor_kernel)
        self.component_counter += 1
        return kappa

    def wall_control(self):
        """ Turns off the filter close to any of the walls in the problem."""
        wall_var = GridVariable('Wall')
        wall_conditions, wall_equations = [], []
        indexes = [OpenSBLIEq(GridVariable('Grid_%d' % direction), self.block.grid_indexes[direction]) for direction in range(self.ndim)]
        wall_equations += indexes
        # Disable the shock filter at any wall boundaries
        for direction in range(self.ndim):
            for side in [0,1]:
                wall = self.wall_boundaries[direction][side]
                if wall:
                    if side == 0:
                        wall_conditions += [ExprCondPair(0, indexes[direction].lhs <= 5)]
                    else:
                        wall_conditions += [ExprCondPair(0, indexes[direction].lhs >= self.block.ranges[direction][side] - 6)]
        # No wall, default condition is the sensor is not turned off
        wall_conditions += [ExprCondPair(1, True)]
        wall_equations += [OpenSBLIEq(wall_var, Piecewise(*wall_conditions))]
        return wall_var, wall_equations

    def filter_application(self, solution_vector, block):
        """ Applies the non-linear filter by subtracting from the q vector after a full RK time-step."""
        resid_kernel = self.residual_kernels[0]
        modified_equations = []
        # Global parameter to control the dissipation to give extra control of the dissipation in the C code
        FC = ConstantObject('shock_filter_control')
        FC.value = 1 # Default condition has no scaling
        CTD.add_constant(FC)
        # The amount of dissipation to apply, using a local flow sensor
        if self.dissipation_sensor == 'Ducros':
            kappa = self.evaluate_Ducros_sensor(block)
        elif self.dissipation_sensor == 'Constant': # No flow sensor for the dissipation control, only a global parameter
            kappa, kappa_evaluation = 1, []
            if self.Mach_correction:
                kappa_Yee, Mach_eqns = self.evaluate_Yee_Mach_sensor(block)
                kappa *= kappa_Yee
        else:
            raise NotImplementedError("Please enter a valid dissipation sensor option: 'Ducros', or 'Constant'.")

        # Take the maximum of the sensor over nearby grid points
        if self.dissipation_sensor == 'Ducros': 
            formula = kappa
            for direction in range(self.ndim):
                for location in [-2, -1, 0, 1, 2]:
                    formula = Max(formula, increment_dataset(kappa, direction, location))

            kappa_max = GridVariable('kappa_max')
            modified_equations += [OpenSBLIEq(kappa_max, formula)]
        else:
            kappa_max = kappa

        # Turn off the sensor at the walls
        wall_detection, wall_equations = self.wall_control()
        modified_equations += wall_equations

        # Apply the filter
        for i, eqn in enumerate(resid_kernel.equations):
            q = self.solution_vector[i]
            modified_equations.append(OpenSBLIEq(q, q + wall_detection*FC*kappa_max*ConstantObject('dt')*eqn.rhs))

        resid_kernel.equations = modified_equations
        residual_kernel = self.create_kernel('Non-linear filter application', resid_kernel.equations, resid_kernel.halo_ranges, block)
        self.component_counter += 1
        self.add_kernel(residual_kernel)
        return

    def zero_work_arrays(self, block):
        resid_kernel = self.residual_kernels[0]
        zeroed_equations = []
        for i in range((self.ndim+2)*self.ndim):
            zeroed_equations += [OpenSBLIEq(block.location_dataset('wk%d' % i), 0.0)]

        zero_kernel = self.create_kernel('Zero the work arrays', zeroed_equations, resid_kernel.halo_ranges, block)
        self.component_counter += 1
        self.add_kernel(zero_kernel)
        return

    def convert_to_datasets(self, block, equations):
        output_equations = []
        for eqn in flatten(equations):
            output_equations += [eqn.convert_to_datasets(block)]
        return output_equations

    def main(self, scheme_order, block):
        # Counter to order the kernels
        self.component_counter = block.blocknumber*10
        # Create the equations for WENO
        eqn = self.create_weno_equations()
        # Convert the equations to datasets on this block
        self.equations = self.convert_to_datasets(block, eqn)
        # Create a WENO scheme
        WS = LLFWeno(scheme_order, formulation='JS', averaging=SimpleAverage([0, 1]), shock_filter=True)
        self.halo_type = set()
        self.halo_type.add(WS.halotype)
        # Start the discretisation and create residual arrays for the equations
        self.Kernels = []
        self.create_residual_arrays(block)
        CR, solution_vector = WS.discretise(self, block)
        # Q vector
        self.solution_vector = flatten(self.time_advance_arrays)
        # Constituent relations evaluations on the Q vector at the end of the full RK time-step
        self.constituent_relations(block)

        # Zero the work arrays
        # self.zero_work_arrays(block)
        # Create the WENO reconstruction kernels
        reconstruction_kernels = []
        for code_gen_order, ker in enumerate(self.reconstruction_kernels):
            halo_ranges = ker.halo_ranges
            reconstruction_kernels.append(self.create_kernel('WENO reconstruction direction %d' % code_gen_order, ker.equations, halo_ranges, block))
            self.component_counter += 1

        self.add_kernel(reconstruction_kernels)
        # Check if there any wall boundary conditions defined on the block.
        self.detect_wall_boundaries()
        # # Create the residual kernel
        self.filter_application(solution_vector, block)
        return
