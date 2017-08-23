from sympy import *
from opensbli.core.opensblifunctions import *
from opensbli.core.grid import GridVariable
from opensbli.core.opensbliequations import *
from opensbli.utilities.helperfunctions import increment_dataset
from opensbli.core.scheme import Scheme
from opensbli.core.weno_opensbli import EigenSystem, Characteristic, LLFCharacteristic
from .kernel import ConstantsToDeclare as CTD

class TenoHalos(object):
    def __init__(self, order, reconstruction=None):
        if not reconstruction:
            if order == 8:
                k = 4
            else:
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
        elif self.order == 6:
            widths = [3, 3, 3, 4]
        elif self.order == 8:
            widths = [3, 3, 3, 4, 4, 5]
        for i, width in enumerate(widths):
            stencils += [TenoStencil(side, width, i)]
        return stencils
    def generate_func_points(self):
        fn_points = []
        side, order = self.side, self.order
        if side == 1:
            fn_points = [[-1, 0, 1], [0, 1, 2], [-2, -1, 0], [0, 1, 2, 3], [-3, -2, -1, 0], [0, 1, 2, 3, 4]]
        elif side == -1:
            fn_points = [[0, 1, 2], [-1, 0, 1], [1, 2, 3], [-2, -1, 0, 1], [0, 1, 2, 3], [-3, -2, -1, 0, 1]]
        return [fn_points[i] for i in range(order-2)]

    def generate_eno_coefficients(self):
        side = self.side
        if side == 1:
            coeffs = [[Rational(-1,6), Rational(5,6), Rational(2,6)], [Rational(2,6), Rational(5,6), Rational(-1,6)], \
            [Rational(2,6), Rational(-7,6), Rational(11,6)]]
            if self.order in [6,8]:
                coeffs += [[Rational(3,12), Rational(13,12), Rational(-5,12), Rational(1,12)]]
            if self.order == 8:
                coeffs += [[Rational(-3,12), Rational(13,12), Rational(-23,12), Rational(25,12)], \
                            [Rational(12,60), Rational(77,60), Rational(-43,60), Rational(17,60), Rational(-3,60)]]
        elif side == -1:
            coeffs = [[Rational(2,6), Rational(5,6), Rational(-1,6)], [Rational(-1,6), Rational(5,6), Rational(2,6)], \
            [Rational(11,6), Rational(-7,6), Rational(1,3)]]
            if self.order in [6,8]:
                coeffs += [[Rational(1,12), Rational(-5,12), Rational(13,12), Rational(3,12)]]
            if self.order == 8:
                coeffs += [[Rational(3,12), Rational(13,12), Rational(-5,12), Rational(1,12)], \
                            [Rational(-1,20), Rational(17,60), Rational(-43,60), Rational(77,60), Rational(12,60)]]
        return coeffs
    def generate_optimal_coefficients(self):
        # Optimal weights are not reversed for TENO when switching between upwind/downwind biasing
        order = self.order
        if self.optimized:
            if order == 5:
                opt_coeffs = [Rational(11,20), Rational(4,10), Rational(1,20)]
            elif order == 6:
                opt_coeffs = [Rational(231,500), Rational(3,10), Rational(27,500), Rational(23,125)]
            elif order == 8:
                opt_coeffs = [0.4336570089737348, 0.2193140179474722, 0.07144766367542149,0.1302093452983125,\
                              0.03089532735084351, 0.1144766367542177]
        elif not self.optimized:
            if order == 5:
                opt_coeffs = [Rational(6,10), Rational(3,10), Rational(1,10)]
            elif order == 6:
                opt_coeffs = [Rational(9,20), Rational(6,20), Rational(1,20), Rational(4,20)]
            elif order == 8: 
                opt_coeffs = [Rational(30,70), Rational(18, 70), Rational(4,70), Rational(12,70), Rational(1,70), Rational(5,70)]
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
        if self.side == 1:
            # Stencils [0, 1, 2] for TENO5
            ## Factored versions: 0, 1, 2, 3
            smoothness_indicators += [Rational(1,4)*(fns[-1] - fns[1])**2 + Rational(13,12)*(fns[-1]-2*fns[0]+fns[1])**2]
            smoothness_indicators += [Rational(1,4)*(3*fns[0]-4*fns[1]+fns[2])**2 + Rational(13,12)*(fns[0]-2*fns[1]+fns[2])**2]
            smoothness_indicators += [Rational(1,4)*(fns[-2]-4*fns[-1]+3*fns[0])**2 + Rational(13,12)*(fns[-2]-2*fns[-1]+fns[0])**2]
            if self.order in [6,8]:
                smoothness_indicators += [Rational(1,36)*(-11*fns[0]+18*fns[1]-9*fns[2]+2*fns[3])**2 + Rational(13,12)*(2*fns[0]-5*fns[1]+4*fns[2]-fns[3])**2 +\
                                      Rational(781,720)*(-fns[0]+3*fns[1]-3*fns[2]+fns[3])]
            # Add Stencils [4, 5] for TENO8
            if self.order == 8:
                smoothness_indicators += [Rational(1,36)*(-2*fns[-3]+9*fns[-2]-18*fns[-1]+11*fns[0])**2 + \
                                        Rational(13,12)*(-1*fns[-3]+4*fns[-2]-5*fns[-1]+2*fns[0])**2 + \
                                        Rational(781,720)*(-1*fns[-3]+3*fns[-2]-3*fns[-1]+fns[0])**2]
                # Use horner on this one, reduces number of operations
                smoothness_indicators += [horner(Rational(1,144)*(-25*fns[0]+48*fns[1]-36*fns[2]+16*fns[3]-3*fns[4])**2 +\
                                        Rational(13,1728)*(35*fns[0]-104*fns[1]+114*fns[2]-56*fns[3]+11*fns[4])**2 +\
                                        Rational(781,2880)*(-5*fns[0]+18*fns[1]-24*fns[2]+14*fns[3]-3*fns[4])**2 -\
                                        Rational(1,4320)*(35*fns[0]-104*fns[1]+114*fns[2]-56*fns[3]+11*fns[4])*\
                                        (fns[0]-4*fns[1]+6*fns[2]-4*fns[3]+fns[4]) + Rational(32803,30240)*\
                                        (fns[0]-4*fns[1]+6*fns[2]-4*fns[3]+fns[4])**2)]
        elif self.side == -1:
            # Stencils [0, 1, 2] for TENO5
            # Factored versions: 0, 1, 2, 3 [0, 1, 2], [-1, 0, 1], [1, 2, 3], [-2, -1, 0, 1]
            smoothness_indicators += [Rational(1,4)*(fns[0]-fns[2])**2 + Rational(13,12)*(fns[0]-2*fns[1]+fns[2])**2]
            smoothness_indicators += [Rational(1,4)*(fns[-1]-4*fns[0]+3*fns[1])**2 + Rational(13,12)*(fns[-1]-2*fns[0]+fns[1])**2]
            smoothness_indicators += [Rational(1,4)*(3*fns[1]-4*fns[2]+fns[3])**2 + Rational(13,12)*(fns[1]-2*fns[2]+fns[3])**2]
            if self.order in [6,8]:
                smoothness_indicators += [Rational(1,36)*(-2*fns[-2]+9*fns[-1]-18*fns[0]+11*fns[1])**2 + \
                                        Rational(13,12)*(-1*fns[-2]+4*fns[-1]-5*fns[0]+2*fns[1])**2 + \
                                        Rational(781,720)*(-1*fns[-2]+3*fns[-1]-3*fns[0]+fns[1])**2]
            if self.order == 8:
                # [0, 1, 2, 3]
                # smoothness_indicators += [Rational(1,36)*(-2*fns[0]-3*fns[1]+6*fns[2]-fns[3])**2 + Rational(13,12)*(fns[0]-2*fns[1]+fns[2])**2 +\
                #                          + Rational(1043,960)*(-fns[0]+3*fns[1]-3*fns[2]+fns[3])**2 + Rational(1,432)*(-2*fns[0]-3*fns[1]+6*fns[2]-fns[3])*\
                #                          (-fns[0]+3*fns[1]-3*fns[2]+fns[3])]
                smoothness_indicators += [horner(Rational(547,240)*fns[0]**2 - Rational(1261,120)*fns[0]*fns[1]+Rational(961,120)*fns[0]*fns[2] - Rational(247,120)*fns[0]*fns[3]\
                                         +Rational(3443,240)*fns[1]**2 -Rational(2983,120)*fns[1]*fns[2] + Rational(267,40)*fns[1]*fns[3]+Rational(2843,240)*fns[2]**2 \
                                         -Rational(821,120)*fns[2]*fns[3] + Rational(89,80)*fns[3]**2)]
                # Factored , [-3, -2, -1, 0, 1]
                smoothness_indicators += [horner(Rational(11329,2520)*fns[-3]**2-Rational(208501,5040)*fns[-3]*fns[-2]+Rational(121621,1680)*fns[-3]*fns[-1] \
                                          -Rational(288007,5040)*fns[-3]*fns[0]+Rational(86329,5040)*fns[-3]*fns[1]+Rational(482963,5040)*fns[-2]**2 \
                                          -Rational(142033,420)*fns[-2]*fns[-1]+Rational(679229,2520)*fns[-2]*fns[0]-Rational(411487,5040)*fns[-2]*fns[1] \
                                          +Rational(507131,1680)*fns[-1]**2 -Rational(68391,140)*fns[-1]*fns[0]+Rational(252941,1680)*fns[-1]*fns[1] \
                                          +Rational(1020563,5040)*fns[0]**2 -Rational(649501,5040)*fns[0]*fns[1]+Rational(53959,2520)*fns[1]**2)]
        fns_dictionary = fns
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
        self.eps = ConstantObject('eps')
        # Parameter to control the spectral properties of the TENO scheme
        self.CT = ConstantObject('TENO_CT')
        CTD.add_constant(self.CT)
        CTD.add_constant(self.eps)
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
            RV.alpha_evaluated.append((C + (tau_5/(self.eps + RV.smoothness_symbols[r])))**q)
        return

    def create_cutoff_equations(self, RV, TC):
        eqns = []
        kronecker_deltas = [Symbol('delta_%d' % r) for r in range(TC.n_stencils)]
        inv_alpha_sum_symbols = [Symbol('inv_alpha_sum')]
        inv_alpha_sum_evaluated = [S.One/sum(RV.alpha_symbols)]
        for r in range(RV.n_stencils):
            eqns += [ceiling((Max(RV.alpha_symbols[r]*inv_alpha_sum_symbols[0]-self.CT, S.Zero)))]
        RV.inv_alpha_sum_symbols = inv_alpha_sum_symbols
        RV.inv_alpha_sum_evaluated = inv_alpha_sum_evaluated
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
        self.eps = ConstantObject('eps')
        # Parameter to control the spectral properties of the TENO scheme
        self.CT = ConstantObject('TENO_CT')
        CTD.add_constant(self.CT)
        CTD.add_constant(self.eps)
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
            RV.alpha_evaluated.append((C + (tau_6/(self.eps + RV.smoothness_symbols[r])))**q)
        return

    def create_cutoff_equations(self, RV, TC):
        eqns = []
        kronecker_deltas = [Symbol('delta_%d' % r) for r in range(TC.n_stencils)]
        inv_alpha_sum_symbols = [Symbol('inv_alpha_sum')]
        inv_alpha_sum_evaluated = [S.One/sum(RV.alpha_symbols)]
        for r in range(RV.n_stencils):
            eqns += [ceiling((Max(RV.alpha_symbols[r]*inv_alpha_sum_symbols[0]-self.CT, S.Zero)))]
        RV.inv_alpha_sum_symbols = inv_alpha_sum_symbols
        RV.inv_alpha_sum_evaluated = inv_alpha_sum_evaluated
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

class Teno8(object):
    def __init__(self):
        # Epsilon to avoid division by zero in non-linear weights
        self.eps = ConstantObject('eps')
        # Parameter to control the spectral properties of the TENO scheme
        self.CT = ConstantObject('TENO_CT')
        CTD.add_constant(self.CT)
        CTD.add_constant(self.eps)
        return

    def __hash__(self):
        h = hash(self._hashable_content())
        self._mhash = h
        return h
    def _hashable_content(self):
        return str(type(cls).__name__)

    def global_smoothness_indicator(self, RV):
        # Global smoothness indicator used in tau_8 for TENO8
        points = [-3, -2, -1, 0, 1, 2, 3, 4]
        f = IndexedBase('f')
        symbolic_functions = []
        fns = {}
        for p in points:
            symbolic_functions.append(f[p])
            fns[p] = symbolic_functions[-1]
        tau_8 = fns[4]*(75349098471*fns[4]-1078504915264*fns[3]+3263178215782*fns[2]\
                -5401061230160*fns[1]+5274436892970*fns[0]-3038037798592*fns[-1]+956371298594*fns[-2]-127080660272*fns[-3])+fns[3]*(3944861897609*fns[3]-24347015748304*fns[2]\
                +41008808432890*fns[1]-40666174667520*fns[0]+23740865961334*fns[-1]-7563868580208*fns[-2]\
                +1016165721854*fns[-3])+fns[2]*(38329064547231*fns[2]-131672853704480*fns[1]+132979856899250*fns[0]-78915800051952*fns[-1]+25505661974314*fns[-2]- 3471156679072*fns[-3])\
                +fns[1]*(115451981835025*fns[1]-238079153652400*fns[0]+144094750348910*fns[-1] \
                    -47407534412640*fns[-2]+6553080547830*fns[-3])+fns[0]*(125494539510175*fns[0] \
                    -155373333547520*fns[-1]+52241614797670*fns[-2]-7366325742800*fns[-3]) \
                    +fns[-1]*(49287325751121*fns[-1]-33999931981264*fns[-2]+4916835566842*fns[-3]) \
                    +fns[-2]*(6033767706599*fns[-2]-1799848509664*fns[-3])+139164877641*fns[-3]*fns[-3]
        tau_8 = Rational(1,62270208000)*(tau_8)
        RV.tau_8_symbol = [Symbol('tau_8')]
        RV.tau_8_evaluated = [tau_8]
        return
    
    def generate_alphas(self, RV, TC):
        """ Create the alpha terms for the non-linear TENO weights.
        arg: object RV: The reconstruction variable object.
        arg: object TenoConfig: Configuration settings for a reconstruction of either left or right."""
        # Scale separation parameters 
        C, q = S.One, 6
        RV.tau_8_symbol = [Symbol('tau')]
        for r in range(TC.n_stencils):
            RV.alpha_symbols += [Symbol('alpha_%d' % r)]
            RV.alpha_evaluated.append((C + (Abs(RV.tau_8_symbol[0])/(self.eps + RV.smoothness_symbols[r])))**q)
        return

    def create_cutoff_equations(self, RV, TC):
        eqns = []
        kronecker_deltas = [Symbol('delta_%d' % r) for r in range(TC.n_stencils)]
        inv_alpha_sum_symbols = [Symbol('inv_alpha_sum')]
        inv_alpha_sum_evaluated = [S.One/sum(RV.alpha_symbols)]
        for r in range(RV.n_stencils):
            eqns += [ceiling((Max(RV.alpha_symbols[r]*inv_alpha_sum_symbols[0]-self.CT, S.Zero)))]
        RV.inv_alpha_sum_symbols = inv_alpha_sum_symbols
        RV.inv_alpha_sum_evaluated = inv_alpha_sum_evaluated
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
        self.inv_alpha_sum_symbols = []
        self.inv_alpha_sum_evaluated = []
        self.tau_8_symbol = []
        self.tau_8_evaluated = []
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
        self.inv_alpha_sum_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.inv_alpha_sum_symbols]
        self.omega_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.omega_symbols]
        self.kronecker_symbols += [GridVariable('%s_%s' % (s, self.name)) for s in original.kronecker_symbols]

        if original.order == 8:
            new_name = self.name.split('_')[-1]
            self.tau_8_symbol = [GridVariable('%s_%s' %(original.tau_8_symbol[0], new_name))]
            originals = original.tau_8_symbol
            new = self.tau_8_symbol[:]
        else:
            originals = []
            new = []
        originals += original.smoothness_symbols + original.alpha_symbols+ original.inv_alpha_sum_symbols + original.kronecker_symbols +original.omega_symbols
        new += self.smoothness_symbols  + self.alpha_symbols+ self.inv_alpha_sum_symbols  + self.kronecker_symbols + self.omega_symbols
        subs_dict = dict(zip(originals, new))
        fn_subs_dict = {}

        for key, value in original.function_stencil_dictionary.iteritems():
            subs_dict[value] = self.function_stencil_dictionary[key]

        self.smoothness_indicators = [s.subs(subs_dict) for s in original.smoothness_indicators]
        # Only update tau_8 for right reconstructions
        if original.order == 8:
            if isinstance(self, LeftTenoReconstructionVariable):
                self.tau_8_evaluated = self.tau_8_symbol
            else:
                self.tau_8_evaluated = [s.subs(subs_dict) for s in original.tau_8_evaluated]
        self.alpha_evaluated = [s.subs(subs_dict) for s in original.alpha_evaluated]
        self.inv_alpha_sum_evaluated = [s.subs(subs_dict) for s in original.inv_alpha_sum_evaluated]
        self.omega_evaluated = [s.subs(subs_dict) for s in original.omega_evaluated]
        self.kronecker_evaluated = [s.subs(subs_dict) for s in original.kronecker_evaluated]
        self.reconstructed_expression = original.reconstructed_expression.subs(subs_dict)
        return

    def add_evaluations_to_kernel(self, kernel):
        all_symbols = self.smoothness_symbols + self.tau_8_symbol + self.alpha_symbols + self.inv_alpha_sum_symbols + self.kronecker_symbols + self.omega_symbols 
        all_evaluations = self.smoothness_indicators + self.tau_8_evaluated + self.alpha_evaluated + self.inv_alpha_sum_evaluated  + self.kronecker_evaluated + self.omega_evaluated

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
    """ Main TENO class."""
    def __init__(self, order, **kwargs):
        """ :arg: int order: Numerical order of the TENO scheme (3,5,...). """
        Scheme.__init__(self, "TenoDerivative", order)
        self.schemetype = "Spatial"
        self.order = order
        if order == 5:
            WT = Teno5()
        elif order == 6:
            WT = Teno6()
        elif order == 8:
            WT = Teno8()
        else:
            raise NotImplementedError("Only 5th, 6th and 8th order TENO implemented currently.")
        # Use optimized schemes? 
        optimized = True
        self.halotype = TenoHalos(self.order)
        self.required_constituent_relations_symbols = {}
        # Generate smoothness coefficients and store configurations for left and right reconstructions.
        self.reconstruction_classes = [RightTenoReconstructionVariable('right'), LeftTenoReconstructionVariable('left')]
        # Populate the quantities required by TENO for the left and right reconstruction variable.
        for no, side in enumerate([1, -1]):
            TC = ConfigureTeno(order, side, optimized)
            # pprint(TenoConfig.func_points)
            RV = self.reconstruction_classes[no]
            RV.order = order
            RV.func_points, RV.n_stencils = sorted(set(flatten(TC.fn_points))), TC.n_stencils
            RV.stencil_points, RV.function_stencil_dictionary = TC.unique_fn_points, TC.fn_dictionary
            RV.smoothness_symbols, RV.smoothness_evaluated = TC.smoothness_symbols, TC.smoothness_indicators
            if (self.order == 8):
                WT.global_smoothness_indicator(RV)
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
            pprint(d.reconstructions)
            for rv in d.reconstructions:
                print type(rv)
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
        grouped = {}
        for cd in all_WDS:
            direction = cd.get_direction[0]
            if direction in grouped.keys():
                grouped[direction] += [cd]
            else:
                grouped[direction] = [cd]
        # TODO: check for size of grouped items
        return grouped