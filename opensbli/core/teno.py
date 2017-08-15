from sympy import *
from opensbli.core.opensblifunctions import *
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliequations import *
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.core.scheme import Scheme
from opensbli.core.weno_opensbli import EigenSystem, Characteristic, LLFCharacteristic

def test_optimal_weights():
    """ Testing function, can delete later."""
    sixth_higher_order = Rational(1,60)*(Symbol('u_-2') - 8*Symbol('u_-1') + 37*Symbol('u_0') + 37*Symbol('u_1') - 8*Symbol('u_2') + Symbol('u_3'))
    right_fifth_higher_order = Rational(1,60)*(2*Symbol('u_-2') - 13*Symbol('u_-1') + 47*Symbol('u_0') + 27*Symbol('u_1') - 3*Symbol('u_2'))
    left_fifth_higher_order = Rational(1,60)*(-3*Symbol('u_-1') + 27*Symbol('u_0') + 47*Symbol('u_1') - 13*Symbol('u_2') + 2*Symbol('u_3'))
    order = 5
    CT = ConfigureTeno(order, 1)
    rhs = []
    opt_coeff_symbols = [Symbol('d_%d' % i) for i in range(order-2)]
    for index, stencil in enumerate(CT.stencils):
        symbols = [Symbol('u_%d' % i) for i in stencil.fn_points]
        coefficients = [i for i in stencil.eno_coeffs]
        local_sum = sum([x*y for (x,y) in zip(symbols, coefficients)])
        rhs += [opt_coeff_symbols[index]*local_sum]
    rhs = sum(rhs)
    eqn = rhs - sixth_higher_order
    sixth_opt_values = [Rational(9,20), Rational(6,20), Rational(1,20), Rational(4,20)]
    fifth_opt_values = [Rational(6,10), Rational(3,10), Rational(1,10)]
    # fifth_opt_values = [Rational(1,10), Rational(6,10), Rational(3,10)]
    # opt_values = [Rational(4,20), Rational(1,20), Rational(6,20), Rational(9,20)]
    sub_list = dict([(x,y) for (x,y) in zip(opt_coeff_symbols, sixth_opt_values)])
    # pprint(sub_list)
    # pprint(eqn)
    # print "\n\n\n"
    # pprint(eqn.subs(sub_list))
    return


class TenoHalos(object):
    def __init__(self, order, reconstruction=None):
        if not reconstruction:
            k = 3 ## Halos correct for TENO5, TENO6
            self.halos = [-k, k+1]
        else:
            self.halos = [-1, 1]
        return

    def get_halos(self, side):
        return self.halos[side]

class TenoStencil(object):
    def __init__(self, side, width, stencil_number):
        self.width = width
        self.side = side
        self.stencil_number = stencil_number
        self.fn_points = None
        self.eno_coeffs = None
        self.opt_coeff = None
        self.smoothness_indicator = None
        return

class ConfigureTeno(object):
    def __init__(self, order, side, optimized=False):
        self.order = order
        self.side = side
        self.optimized = optimized
        self.n_stencils = order - 2
        # Create stencils
        self.stencils = self.generate_stencils()
        self.fn_points = self.generate_func_points()
        self.eno_coeffs = self.generate_eno_coefficients()
        self.opt_coeffs = self.generate_optimal_coefficients()
        self.fn_dictionary, self.smoothness_indicators, self.smoothness_symbols = self.generate_smoothness_indicators()
        self.unique_fn_points = sorted(set(flatten(self.fn_points)))
        # Add the function points and coefficients to the stencil objects
        self.update_stencils()
        return 

    def generate_stencils(self):
        stencils = []
        side = self.side
        if self.order == 5:
            widths = [3, 3, 3]
            for i, width in enumerate(widths):
                stencils += [TenoStencil(side, width, i)]
        elif self.order == 6:
            widths = [3, 3, 3, 4]
            for i, width in enumerate(widths):
                stencils += [TenoStencil(side, width, i)]
        return stencils
    def generate_func_points(self):
        fn_points = []
        side = self.side
        order = self.order
        if side == 1:
            fn_points = [[-1, 0, 1], [0, 1, 2], [-2, -1, 0]]
            if order == 6:
                fn_points += [[0, 1, 2, 3]]
        elif side == -1:
            fn_points = [[0, 1, 2], [-1, 0, 1], [1, 2, 3]]
            if order == 6:
                fn_points += [[-2, -1, 0, 1]]
        return fn_points

    def generate_eno_coefficients(self):
        side = self.side
        if side == 1:
            coeffs = [[Rational(-1,6), Rational(5,6), Rational(2,6)], [Rational(2,6), Rational(5,6), Rational(-1,6)], \
            [Rational(2,6), Rational(-7,6), Rational(11,6)]]
            if self.order == 6:
                coeffs += [[Rational(3,12), Rational(13,12), Rational(-5,12), Rational(1,12)]]
        elif side == -1:
            coeffs = [[Rational(2,6), Rational(5,6), Rational(-1,6)], [Rational(-1,6), Rational(5,6), Rational(2,6)], \
            [Rational(11,6), Rational(-7,6), Rational(1,3)]]
            if self.order == 6:
                coeffs += [[Rational(1,12), Rational(-5,12), Rational(13,12), Rational(3,12)]]
        return coeffs
    def generate_optimal_coefficients(self):
        # Optimal weights are not reversed for TENO when switching between upwind/downwind biasing
        order = self.order
        if self.optimized:
            if order == 6:
                opt_coeffs = [Rational(231,500), Rational(3,10), Rational(27,500), Rational(23,125)]
            elif order == 5:
                opt_coeffs = [Rational(11,20), Rational(4,10), Rational(1,20)]

        elif not self.optimized:
            if order == 6:
                opt_coeffs = [Rational(9,20), Rational(6,20), Rational(1,20), Rational(4,20)]
            elif order == 5:
                opt_coeffs = [Rational(6,10), Rational(3,10), Rational(1,10)]
        return opt_coeffs
    def generate_smoothness_indicators(self):
        points = sorted(list(set(flatten([i for i in self.fn_points]))))
        f = IndexedBase('f')
        symbolic_functions = []
        fns = {}
        for p in points:
            symbolic_functions.append(f[p])
            fns[p] = symbolic_functions[-1]
        smoothness_indicators = []
        # 5th order stencils
        if self.side == 1:
            smoothness_indicators += [4*fns[-1]**2 -13*fns[-1]*fns[0] + 5*fns[-1]*fns[1] + 13*fns[0]**2 \
                                        -13*fns[0]*fns[1] + 4*fns[1]**2]
            smoothness_indicators += [10*fns[0]**2 - 31*fns[0]*fns[1] + 11*fns[0]*fns[2] + 25*fns[1]**2 \
                                        -19*fns[1]*fns[2] + 4*fns[2]**2]
            smoothness_indicators += [25*fns[-1]**2 -19*fns[-1]*fns[-2] - 31*fns[-1]*fns[0] + 4*fns[-2]**2 \
                                        +11*fns[-2]*fns[0] + 10*fns[0]**2]
        elif self.side == -1:
            smoothness_indicators += [4*fns[0]**2 - 13*fns[0]*fns[1] + 5*fns[0]*fns[2] + 13*fns[1]**2 \
                                    - 13*fns[1]*fns[2] + 4*fns[2]**2]
            smoothness_indicators += [4*fns[-1]**2 - 19*fns[-1]*fns[0] + 11*fns[-1]*fns[1] + 25*fns[0]**2 \
                                        -31*fns[0]*fns[1] + 10*fns[1]**2]
            smoothness_indicators += [10*fns[1]**2 - 31*fns[1]*fns[2] + 11*fns[1]*fns[3]+25*fns[2]**2 \
                                        -19*fns[2]*fns[3] + 4*fns[3]**2]
        smoothness_indicators = [Rational(1,3)*i for i in smoothness_indicators]
        # Add extra smoothness indicator for 6th order
        if self.order == 6:
            if self.side == 1:
                four_point_indicator = [547*fns[0]**2 - 2522*fns[0]*fns[1] + 1922*fns[0]*fns[2]-494*fns[0]*fns[3] \
                                        +3443*fns[1]**2 - 5966*fns[1]*fns[2] + 1602*fns[1]*fns[3] + 2843*fns[2]**2 \
                                        - 1642*fns[2]*fns[3] + 267*fns[3]**2]
                smoothness_indicators += [Rational(1,240)*i for i in four_point_indicator]
            elif self.side == -1:
                four_point_indicator = [7043*fns[-1]**2 - 3882*fns[-1]*fns[-2] -17246*fns[-1]*fns[0] \
                                        +7042*fns[-1]*fns[1] + 547*fns[-2]**2 + 4642*fns[-2]*fns[0] \
                                        - 1854*fns[-2]*fns[1] + 11003*fns[0]**2 - 9402*fns[0]*fns[1] + 2107*fns[1]**2]
                smoothness_indicators += [Rational(1,240)*i for i in four_point_indicator]
        fns_dictionary = fns
        # for item in smoothness_indicators:
        #   pprint(item)
        #   pprint(item.count_ops())
        #   new_item = horner(item)
        #   pprint(new_item)
        #   pprint(new_item.count_ops())
        smoothness_indicators = [horner(eqn) for eqn in smoothness_indicators]
        smoothness_symbols = [Symbol('beta_%d' % r) for r in range(len(smoothness_indicators))]
        return fns_dictionary, smoothness_indicators, smoothness_symbols

    def update_stencils(self):
        for stencil in self.stencils:
            no = stencil.stencil_number
            stencil.fn_points = self.fn_points[no]
            stencil.eno_coeffs = self.eno_coeffs[no]
            stencil.opt_coeff = self.opt_coeffs[no]
            stencil.smoothness_indicator = self.smoothness_indicators[no]
        return

class Teno5(object):
    def __init__(self):
        # Epsilon to avoid division by zero in non-linear weights
        self.eps = 1e-15
        # Parameter to control the spectral properties of the TENO scheme
        self.CT = 1e-5
        return
    
    def generate_alphas(self, RV, TC):
        """ Create the alpha terms for the non-linear TENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        # Scale separation parameters 
        C, q = S.One, 6
        # Global reference smoothness indicator tau_5
        tau_5 = Abs(RV.smoothness_symbols[0] - RV.smoothness_symbols[-1])
        for r in range(TC.n_stencils):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append((C + (tau_5/(self.eps + RV.smoothness_symbols[r]))**q))
        return
    def create_cutoff_equations(self, RV, TC):
        eqns = []
        kronecker_deltas = [Symbol('delta_%d' % r) for r in range(TC.n_stencils)]
        chi_symbols = [Symbol('Chi_%d' % r) for r in range(TC.n_stencils)]
        chi_evaluated = [RV.alpha_symbols[r]/sum(RV.alpha_symbols) for r in range(TC.n_stencils)]

        for r in range(RV.n_stencils):
            ecs = [ExprCondPair(S.Zero, StrictLessThan(chi_symbols[r],self.CT))]
            ecs += [ExprCondPair(S.One, True)]
            eqns += [Piecewise(*ecs, **{'evaluate':False})]
        RV.chi_symbols = chi_symbols
        RV.chi_evaluated = chi_evaluated
        RV.kronecker_evaluated = eqns
        RV.kronecker_symbols = kronecker_deltas
        return

    def generate_omegas(self, RV, TC):
        """ Create the omega terms for the non-linear TENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        weight_sum = sum([TC.opt_coeffs[i]*RV.kronecker_symbols[i] for i in range(RV.n_stencils)])
        RV.omega_symbols = [Symbol('omega_%d' % r) for r in range(RV.n_stencils)]
        RV.omega_evaluated = [TC.opt_coeffs[r]*RV.kronecker_symbols[r] / weight_sum for r in range(RV.n_stencils)]
        return

    def generate_reconstruction(self, RV, TenoConfig):
        """ Create the final TENO stencil by summing the stencil points, ENO coefficients and TENO weights..
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        reconstruction = 0
        fns = []
        for stencil in TenoConfig.stencils:
            RV.smoothness_indicators.append(stencil.smoothness_indicator)
            fns = [RV.function_stencil_dictionary[i] for i in stencil.fn_points]
            eno_interpolation = sum([point*coefficient for (point, coefficient) in zip(fns, stencil.eno_coeffs)])
            reconstruction += RV.omega_symbols[stencil.stencil_number]*eno_interpolation
        RV.reconstructed_expression = reconstruction
        return

class Teno6(object):
    def __init__(self):
        # Epsilon to avoid division by zero in non-linear weights
        self.eps = 1e-15
        # Parameter to control the spectral properties of the TENO scheme
        self.CT = 1e-5
        return
    
    def generate_alphas(self, RV, TC):
        """ Create the alpha terms for the non-linear TENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        # Scale separation parameters 
        C, q = S.One, 6
        # Global reference smoothness indicator tau_6
        tau_6 = RV.smoothness_symbols[3] - Rational(1,6)*(RV.smoothness_symbols[0] + RV.smoothness_symbols[2] - 4*RV.smoothness_symbols[1])
        for r in range(TC.n_stencils):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append((C + (tau_6/(self.eps + RV.smoothness_symbols[r]))**q))
        return
    def create_cutoff_equations(self, RV, TC):
        eqns = []
        kronecker_deltas = [Symbol('delta_%d' % r) for r in range(TC.n_stencils)]
        chi_symbols = [Symbol('Chi_%d' % r) for r in range(TC.n_stencils)]
        chi_evaluated = [RV.alpha_symbols[r]/sum(RV.alpha_symbols) for r in range(TC.n_stencils)]

        for r in range(RV.n_stencils):
            ecs = [ExprCondPair(0, StrictLessThan(chi_symbols[r],self.CT))]
            ecs += [ExprCondPair(1, True)]
            eqns += [Piecewise(*ecs, **{'evaluate':False})]
        RV.chi_symbols = chi_symbols
        RV.chi_evaluated = chi_evaluated
        RV.kronecker_evaluated = eqns
        RV.kronecker_symbols = kronecker_deltas
        return

    def generate_omegas(self, RV, TC):
        """ Create the omega terms for the non-linear TENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        weight_sum = sum([TC.opt_coeffs[i]*RV.kronecker_symbols[i] for i in range(RV.n_stencils)])
        RV.omega_symbols = [Symbol('omega_%d' % r) for r in range(RV.n_stencils)]
        RV.omega_evaluated = [TC.opt_coeffs[r]*RV.kronecker_symbols[r] / weight_sum for r in range(RV.n_stencils)]
        return

    def generate_reconstruction(self, RV, TenoConfig):
        """ Create the final TENO stencil by summing the stencil points, ENO coefficients and TENO weights..
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        reconstruction = 0
        fns = []
        for stencil in TenoConfig.stencils:
            RV.smoothness_indicators.append(stencil.smoothness_indicator)
            fns = [RV.function_stencil_dictionary[i] for i in stencil.fn_points]
            eno_interpolation = sum([point*coefficient for (point, coefficient) in zip(fns, stencil.eno_coeffs)])
            reconstruction += RV.omega_symbols[stencil.stencil_number]*eno_interpolation
        RV.reconstructed_expression = reconstruction
        return


class TenoReconstructionVariable(object):
    def __init__(self, name):
        """ Reconstruction variable object to hold the quantities required for TENO.
        arg: str: name: Name of the reconstruction, either left or right."""
        self.name = name
        self.smoothness_indicators = []
        self.smoothness_symbols = []
        self.alpha_evaluated = []
        self.alpha_symbols = []
        self.chi_evaluated = []
        self.chi_symbols = []
        self.omega_evaluated = []
        self.omega_symbols = []
        self.kronecker_evaluated = []
        self.kronecker_symbols = []
        self.function_stencil_dictionary = {}
        self.reconstructed_expression = None
        self.reconstructed_symbol = GridVariable('%s_%s' % ('reconstruct', name))
        return

    def update_quantities(self, original):
        """ Updates the quantities required by TENO in the reconstruction variable.
        arg: object: original: Reconstruction object variable, either left or right reconstruction."""
        print "Reconstruction variable: ", self.name
        self.smoothness_symbols += [GridVariable('%s%s' % (s, self.name)) for s in original.smoothness_symbols]
        self.alpha_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.alpha_symbols]
        self.chi_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.chi_symbols]
        self.omega_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.omega_symbols]
        self.kronecker_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.kronecker_symbols]


        subs_dict = dict(zip(original.smoothness_symbols + original.alpha_symbols+ original.chi_symbols + original.kronecker_symbols +original.omega_symbols, \
                             self.smoothness_symbols + self.alpha_symbols+ self.chi_symbols + self.kronecker_symbols + self.omega_symbols))
        fn_subs_dict = {}
        for key, value in original.function_stencil_dictionary.iteritems():
            print key, value
            subs_dict[value] = self.function_stencil_dictionary[key]
        self.smoothness_indicators = [s.subs(subs_dict) for s in original.smoothness_indicators]
        self.alpha_evaluated = [s.subs(subs_dict) for s in original.alpha_evaluated]
        self.chi_evaluated = [s.subs(subs_dict) for s in original.chi_evaluated]
        self.omega_evaluated = [s.subs(subs_dict) for s in original.omega_evaluated]
        self.kronecker_evaluated = [s.subs(subs_dict) for s in original.kronecker_evaluated]
        self.reconstructed_expression = original.reconstructed_expression.subs(subs_dict)
        return

    def add_evaluations_to_kernel(self, kernel):
        all_symbols = self.smoothness_symbols + self.alpha_symbols + self.chi_symbols + self.kronecker_symbols + self.omega_symbols 
        all_evaluations = self.smoothness_indicators + self.alpha_evaluated + self.chi_evaluated + self.kronecker_evaluated + self.omega_evaluated
        for no, value in enumerate(all_symbols):
            kernel.add_equation(Eq(value, all_evaluations[no]))
        kernel.add_equation(Eq(self.reconstructed_symbol, self.reconstructed_expression))
        return


class LeftTenoReconstructionVariable(TenoReconstructionVariable):
    def __init__(self, name):
        """ Reconstruction object for the left reconstruction.
        arg: str: name: 'left' """
        TenoReconstructionVariable.__init__(self, name)
        return


class RightTenoReconstructionVariable(TenoReconstructionVariable):
    def __init__(self, name):
        """ Reconstruction object for the right reconstruction.
        arg: str: name: 'right' """
        TenoReconstructionVariable.__init__(self, name)
        return

class Teno(Scheme):
    """ Main TENO class. Performs the Jiang-Shu TENO reconstruction procedure. Refer to the reference:
    'Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory schemes for Hyperbolic Conservation
    laws. Shu (1997). The alpha & omega quantities follow the description from this paper."""
    def __init__(self, order, **kwargs):
        """ :arg: int order: Numerical order of the TENO scheme (3,5,...). """
        Scheme.__init__(self, "TenoDerivative", order)
        self.schemetype = "Spatial"
        if order == 5:
            self.order = 5
            WT = Teno5()
        elif order == 6:
            self.order = 6
            WT = Teno6()
        else:
            raise NotImplementedError("Only 5th and 6th order TENO implemented currently.")
        # Use optimized schemes? 
        optimized = True
        self.halotype = TenoHalos(self.order)
        self.required_constituent_relations_symbols = {}
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        self.reconstruction_classes = [LeftTenoReconstructionVariable('left'), RightTenoReconstructionVariable('right')]
        # Populate the quantities required by TENO for the left and right reconstruction variable.
        for no, side in enumerate([-1, 1]):
            TC = ConfigureTeno(order, side, optimized)
            # pprint(TenoConfig.func_points)
            RV = self.reconstruction_classes[no]
            RV.func_points, RV.n_stencils = sorted(set(flatten(TC.fn_points))), TC.n_stencils
            RV.stencil_points, RV.function_stencil_dictionary = TC.unique_fn_points, TC.fn_dictionary
            RV.smoothness_symbols, RV.smoothness_evaluated = TC.smoothness_symbols, TC.smoothness_indicators
            WT.generate_alphas(RV, TC)
            WT.create_cutoff_equations(RV, TC)
            WT.generate_omegas(RV, TC)
            WT.generate_reconstruction(RV, TC)
            self.reconstruction_classes[no] = RV
        return

    def interpolate_reconstruction_variables(self, derivatives, kernel):
        """ Perform the TENO interpolation on the reconstruction variables.
        arg: list: derivatives: A list of the TENO derivatives to be computed.
        arg: object: kernel: The current computational kernel.
        """
        for d in derivatives:
            pprint(d)
            for rv in d.reconstructions:
                if isinstance(rv, type(self.reconstruction_classes[1])):
                    original_rv = self.reconstruction_classes[1]
                elif isinstance(rv, type(self.reconstruction_classes[0])):
                    original_rv = self.reconstruction_classes[0]
                else:
                    raise ValueError("Reconstruction must be left or right")
                rv.update_quantities(original_rv)
                rv.add_evaluations_to_kernel(kernel)
        return

    def update_constituent_relation_symbols(self, sym, direction):
        """ Function to take the set of required quantities from the constituent relations in symbolic form
         and update the directions in which they are used.
        arg: set: sym: Set of required symbols.
        arg: int: direction: The axis on which TENO is being applied to (x0, x1 ..)."""
        if isinstance(sym, Symbol):
            sym = [sym]
        elif isinstance(sym, list) or isinstance(sym, set):
            pass
        else:
            raise ValueError("The symbol provided should be either a list or symbols")

        for s in sym:
            if isinstance(s, DataSetBase):
                s = s.noblockname
            if s in self.required_constituent_relations_symbols.keys():
                self.required_constituent_relations_symbols[s] += [direction]
            else:
                self.required_constituent_relations_symbols[s] = [direction]
        return

    def generate_constituent_relations_kernels(self, block):
        """ Generates constituent relation kernels on the block.
        arg: object: block: The current block."""
        crs = {}
        for key in self.required_constituent_relations_symbols:
            pprint([key, self.required_constituent_relations_symbols[key]])
            kernel = Kernel(block, computation_name="CR%s" % key)
            kernel.set_grid_range(block)
            for direction in self.required_constituent_relations_symbols[key]:
                kernel.set_halo_range(direction, 0, self.halotype)
                kernel.set_halo_range(direction, 1, self.halotype)
            crs[DataSetBase(str(key))[DataSetBase.location()]] = kernel
        return crs

class LLFTeno(LLFCharacteristic, Teno):
    def __init__(self, eigenvalue, left_ev, right_ev, order, ndim, averaging=None):
        LLFCharacteristic.__init__(self, eigenvalue, left_ev, right_ev, order, ndim, averaging)
        Teno.__init__(self, order)
        return

    def reconstruction_halotype(self, order, reconstruction=True):
        return TenoHalos(order,reconstruction)
    def group_by_direction(self, eqs):
        """ Groups the input equations by the direction (x0, x1, ...) they depend upon.
        arg: list: eqs: List of equations to group by direction.
        returns: dict: grouped: Dictionary of {direction: equations} key, value pairs for equations grouped by direction."""
        all_WDS = []
        for eq in eqs:
            all_WDS += list(eq.atoms(TenoDerivative))
            # pprint([eq, eq.atoms(TenoDerivative)])
        # all_WDS = list(set(all_WDS))
        # pprint(all_WDS)
        # if len(all_WDS) != len(eqs):
        #     raise ValueError("Number of WD in each equation should be one.")
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        # TODO: check for size of grouped items
        return grouped