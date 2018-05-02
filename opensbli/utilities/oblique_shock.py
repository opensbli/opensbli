from sympy import Symbol, sin, Eq, tan, sqrt, atan, cot, cos, pi, pprint, rad, asin


class ObliqueShock(object):
    """ Class to compute the oblique shock relations.

    :arg float wave_angle: Oblique shock angle beta in degrees.
    :arg float input_mach_number: Free-stream pre-shock condition Mach number.
    :arg float gamma: Ratio of gas constants."""

    def __init__(self, wave_angle, input_mach_number, gamma):
        print "Input Mach number and wave angle are %f, %f \n" % (input_mach_number, wave_angle)
        # Check that wave angle is greater than the Mach angle
        if (wave_angle < asin(1.0/input_mach_number)):
            raise ValueError("Wave angle must be greater than the Mach angle.")
        self.M1, self.beta, self.gamma = Symbol('M_1'), Symbol('beta'), Symbol('gamma')
        symbols = [self.M1, self.beta, self.gamma]
        values = [input_mach_number, rad(wave_angle), gamma]
        # Substitution dictionary for the input symbols/values
        self.subs_dict = dict([(x, y) for (x, y) in zip(symbols, values)])
        wedge_angle = self.theta_beta_M_eqn()
        p2p1 = self.pressure_ratio()
        rho2rho1 = self.density_ratio()
        T2T1 = self.temperature_ratio(rho2rho1, p2p1)
        M2 = self.post_shock_mach_number(wedge_angle)
        symbols = [Symbol('wave_angle'), Symbol('wedge_angle'), Symbol('rho2rho1'), Symbol('p2p1'), Symbol('T2T1'), Symbol('M_2')]
        values = [rad(wave_angle).evalf(), wedge_angle, rho2rho1, p2p1, T2T1, M2]
        self.subs_dict.update(dict([(x, y) for (x, y) in zip(symbols, values)]))
        return

    def theta_beta_M_eqn(self):
        """ Function to calculate the wedge angle for this Mach number and wave angle. 
        Evaluates the equation tan(theta) = 2*cot(beta)*(M1**2 * sin^2(beta) - 1)/(M1**2*(gama + cos(2beta)) + 2). """
        beta, M1, gamma = self.beta, self.M1, self.gamma
        tbm_eqn = Eq(Symbol('theta'), atan(2*cot(beta)*(M1**2 * sin(beta)**2 - 1)/(M1**2 * (gamma+cos(2*beta))+2)))
        return tbm_eqn.subs(self.subs_dict).evalf().rhs

    def pressure_ratio(self):
        """ Evaluates pressure ratio p2/p1 for a given input Mach number and wave angle."""
        beta, M1, gamma = self.beta, self.M1, self.gamma
        p2_p1 = Eq(Symbol('p2p1'), 1 + (2*gamma)/(gamma+1) * (M1**2 * sin(beta)**2 - 1))
        return p2_p1.subs(self.subs_dict).evalf().rhs

    def density_ratio(self):
        """ Evaluates density ratio rho2/rho1 for a given input Mach number and wave angle."""
        beta, M1, gamma = self.beta, self.M1, self.gamma
        rho2_rho1 = Eq(Symbol('rho2rho1'), ((gamma+1)*(M1**2 * sin(beta)**2))/((gamma-1)*M1**2 * sin(beta)**2 + 2))
        return rho2_rho1.subs(self.subs_dict).evalf().rhs

    def temperature_ratio(self, density_ratio, pressure_ratio):
        """ Evaluates temperature ratio T2/T1 for a given input Mach number and wave angle."""
        T2_T1 = Eq(Symbol('T2T1'), pressure_ratio*(1.0/(density_ratio)))
        return T2_T1.rhs

    def post_shock_mach_number(self, wedge_angle):
        """ Computes the post-shock Mach number M2."""
        beta, M1, gamma, theta = self.beta, self.M1, self.gamma, wedge_angle
        M2 = Symbol('M_2')
        M2_eqn = Eq(M2, (1/sin(beta-theta))*(sqrt((1+0.5*(gamma-1)*M1**2 * sin(beta)**2)/(gamma*M1**2 * sin(beta)**2 - 0.5*(gamma-1)))))
        return M2_eqn.subs(self.subs_dict).evalf().rhs


class ShockConditions(ObliqueShock):
    """ Class to calculate oblique shock relations for conservative/primitive variables.

    :arg float wave_angle: Oblique shock angle beta in degrees.
    :arg float input_mach_number: Free-stream pre-shock condition Mach number.
    :arg float gamma: Ratio of gas constants."""
    def __init__(self, wave_angle, input_mach_number, gamma, pre_shock_conditions=None):
        ObliqueShock.__init__(self, wave_angle, input_mach_number, gamma)
        # TO DO: Finish testing pre_shock_conditions part, and solve for theta using secant method
        return

    def conservative_post_shock_conditions(self, freestream_velocity):
        wave, wedge = Symbol('wave_angle'), Symbol('wedge_angle')
        rho2rho1, p2p1, T2T1 = Symbol('rho2rho1'), Symbol('p2p1'), Symbol('T2T1')
        M1, M2, gamma = Symbol('M_1'), Symbol('M_2'), Symbol('gamma')
        free_stream = freestream_velocity
        u = (free_stream*cos(0.5*pi - wave)*tan(wave - wedge)/tan(wave)).subs(self.subs_dict).evalf()
        v = (free_stream*cos(0.5*pi - wave)/tan(wave)).subs(self.subs_dict).evalf()
        rhou = (rho2rho1*sqrt(u*u+v*v)*cos(wedge)).subs(self.subs_dict)
        rhov = (-rho2rho1*sqrt(u*u+v*v)*sin(wedge)).subs(self.subs_dict)
        # Calculate post-shock pressure
        # Freestream pressure for this normalisation.
        p1 = 1.0/(gamma*M1**2)
        p2 = (p2p1*p1).subs(self.subs_dict)
        # Calculate post-shock total energy
        rho = rho2rho1.subs(self.subs_dict)
        rhoE = ((p2/(gamma-1) + (0.5/rho2rho1)*(rhou**2 + rhov**2))).subs(self.subs_dict)
        print "conservative values post shock: rho, rhou, rhov, rhoE"
        print rho, rhou, rhov, rhoE
        return [rho, rhou, rhov, rhoE]

    def primitive_post_shock_conditions(self, freestream_velocity):
        wave, wedge = Symbol('wave_angle'), Symbol('wedge_angle')
        rho2rho1, p2p1, T2T1 = Symbol('rho2rho1'), Symbol('p2p1'), Symbol('T2T1')
        M1, M2, gamma = Symbol('M_1'), Symbol('M_2'), Symbol('gamma')
        free_stream = freestream_velocity
        u = (free_stream*cos(0.5*pi - wave)*tan(wave - wedge)/tan(wave)).subs(self.subs_dict).evalf()
        v = (free_stream*cos(0.5*pi - wave)/tan(wave)).subs(self.subs_dict).evalf()
        u_out= (sqrt(u*u+v*v)*cos(wedge)).subs(self.subs_dict)
        v_out = (-1.0*sqrt(u*u+v*v)*sin(wedge)).subs(self.subs_dict)
        # Calculate post-shock pressure
        # Freestream pressure for this normalisation.
        p1 = 1.0/(gamma*M1**2)
        p2 = (p2p1*p1).subs(self.subs_dict)
        # Calculate post-shock total energy
        rho = rho2rho1.subs(self.subs_dict)
        print "primitive values post shock: rho, u, v, p"
        print rho, u_out, v_out, p2
        return [rho, u_out, v_out, p2]
