from opensbli.utilities.numerical_functions import spline, splint
from sympy import Eq, Piecewise
from scipy.integrate import odeint
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from opensbli.initialisation import GridBasedInitialisation
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.core.grid import GridVariable
import warnings

plt.style.use('classic')


class Boundary_layer_profile(object):
    def __init__(self, xmach, Pr, gama, Tw, Re, length, npoints, beta):
        """ Performs a similarity solution (Viscous fluid flow, F.White 1974),
        to obtain u and T profiles for a laminar compressible boundary-layer.
        arg: float: xmach: Free-stream Mach numnber.
        arg: float: Pr: Prandtl number.
        arg: float: gama: Ratio of specific heats.
        arg: float: Tw: Wall temperature, use -1 for adibatic wall conditions.
        arg: float: Re: Free-stream Reynolds number.
        arg: float: length: Doman length in the wall-normal direction.
        arg: int: npoints: Number of points in the wall-normal direction.
        arg: float: beta: Value of the stretching factor for non-uniform grids.
        """
        self.y, self.u, self.T, self.scale = self.generate_boundary_layer_profile(xmach, Pr, gama, Tw, Re)
        self.Re = Re
        self.interpolate(length, npoints, beta)
        self.n = np.size(self.y)
        return

    def interpolate(self, domain_length, npoints, beta):
        """ Interpolates the boundary-layer data to new coordinates.
        arg: float: domain_length: Length of the new coordinate domain.
        arg: int: npoints: Number of points in the new coordinates.
        arg: float: beta: Value of the stretching factor for stretched coordinates.
        returns: None. """
        eta = np.linspace(0, 1, npoints)
        # y_full = np.linspace(0, domain_length, npoints)
        y_full = domain_length*np.sinh(beta*eta)/np.sinh(beta)
        y, u, T = self.y[:], self.u[:], self.T[:]
        n = self.n
        d2y_u = spline(y, u, n, 0, 0)
        d2y_T = spline(y, T, n, 0, 0)
        u_full, T_full = np.zeros_like(y_full), np.zeros_like(y_full)
        n_full = np.size(u_full)
        for i in range(n_full):
            if (y_full[i] < y[-1]):
                u_full[i] = splint(y, u, d2y_u, n, y_full[i])
                T_full[i] = splint(y, T, d2y_T, n, y_full[i])
            else:  # Stop spline interpolation if we are at the edge of the boundary-layer
                u_full[i] = u_full[i-1]
                T_full[i] = T_full[i-1]
        self.y, self.u, self.T = y_full, u_full, T_full
        return

    def compbl(self, v, p=None):
        """ Sets up the system of equations to be integrated by odeint.
        arg: ndarray: v: Solution vector.
        arg: None: p: Empty dummy argument required by the odeint function.
        returns: list: dv: System of equations.
        """
        suth = self.suth
        c = np.sqrt(v[3])*(1.0+suth)/(v[3]+suth)
        dcdg = 1.0/(2.0*np.sqrt(v[3])) - np.sqrt(v[3])/(v[3]+suth)
        dcdg *= (1.0+suth) / (v[3]+suth)
        cp = dcdg*v[4]
        dv = [v[1], v[2], -v[2]*(cp+v[0])/c, v[4],
              -v[4]*(cp+self.pr*v[0])/c - self.pr*(self.gama-1)*self.xmach**2 * v[2]**2]
        return dv

    def generate_boundary_layer_profile(self, xmach, pr, gama, Tw, Re):
        """ Generates a boundary layer initial profile. Solves the mean flow
        in a compressible boundary layer. (Equations 7.32) in White (1974).
        arg txtfile: Text file containing stretched y-coordinates.
        arg float: xmach: Mach number.
        arg float: pr: Prandtl number.
        arg float: gama: Ratio of specific heats.
        arg float: Tw: Wall temperature Tw/Tinf (< 0 for adiabatic)
        arg float: Re: Reynolds number."""
        n_iter, jmax = 5, 1001
        v, dv, f, f1, f2 = np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)
        self.pr, self.gama, self.xmach, self.Re, self.Tw = pr, gama, xmach, Re, Tw
        self.nvisc = 1  # 1 Sutherlands law, 2 Power law, 3 Chapman-Rubesin approximation.
        Tinf = 288.0
        etamax = 10.0
        nstep = jmax-1
        errtol = 1e-10
        self.suth = 110.40/Tinf  # Sutherland constant
        self.soln = np.zeros((5, jmax))
        self.eta = np.linspace(0, etamax, jmax)
        vstart = np.zeros(5)
        if Tw > 0:
            vstart[3] = Tw
        elif Tw < 0:  # Adiabatic wall, initial guess
            v[1] = 1.0 + 0.5*0.84*(gama-1)*xmach**2
            v[0] = 0.47*v[1]**0.21
        else:  # Fixed wall temperature, initial guess. Tested for 1<Tw<4 and M<8 (Pr=0.72, gamma=1.4)
            v[1] = 4.5  # !0.062*xmach**2-0.1*(Tw-1.0)*(1.0+xmach)/(0.2+xmach)
            v[0] = 0.494891  # !0.45-0.01*xmach+(Tw-1.)*0.06
        # Initial increments
        dv[0] = v[0]*0.01
        if Tw < 0:
            dv[1] = v[1]*0.01
        else:
            dv[1] = 0.1
        # Main loop
        # ODE solver parameters
        abserr = 1.0e-8
        relerr = 1.0e-6
        for k in range(n_iter):
            vstart[2] = v[0]  # Initial value
            if Tw < 0:
                vstart[3] = v[1]
            else:
                vstart[4] = 1.141575  # v[1]
            # Call the ODE solver.
            self.soln = odeint(self.compbl, vstart, self.eta,
                               atol=abserr, rtol=relerr).T
            err = np.abs(self.soln[1, nstep]-1.0) + np.abs(self.soln[3, nstep]-1.0)
            if err < errtol:
                break
            f[0] = self.soln[1, nstep] - 1.0
            f[1] = self.soln[3, nstep] - 1.0
            # Increment v[0]
            vstart[2] = v[0] + dv[0]
            if (Tw < 0):
                vstart[3] = v[1]
            else:
                vstart[4] = v[1]
            self.soln = odeint(self.compbl, vstart, self.eta,
                               atol=abserr, rtol=relerr).T
            f1[0] = self.soln[1, nstep] - 1.0
            f1[1] = self.soln[3, nstep] - 1.0
            # Increment v[1]
            vstart[2] = v[0]
            if (Tw < 0):
                vstart[3] = v[1] + dv[1]
            else:
                vstart[4] = v[1] + dv[1]
            self.soln = odeint(self.compbl, vstart, self.eta,
                               atol=abserr, rtol=relerr).T
            f2[0] = self.soln[1, nstep] - 1.0
            f2[1] = self.soln[3, nstep] - 1.0
            # Solve the linear system
            al11 = (f1[0] - f[0])/dv[0]
            al21 = (f1[1] - f[1])/dv[0]
            al12 = (f2[0] - f[0])/dv[1]
            al22 = (f2[1] - f[1])/dv[1]
            det = float(al11*al22 - al21*al12)
            # New dv for improved solution
            dv[0] = (-al22*f[0] + al12*f[1])/det
            dv[1] = (al21*f[0] - al11*f[1])/det
            v[0] += dv[0]
            v[1] += dv[1]
            # Write out improved estimate
            print 'it   ', k, dv[0], dv[1], v[0], v[1]
        print "final vaues ", v[0], v[1], err
        y, u, T, scale = self.integrate_boundary_layer(nstep)
        return y, u, T, scale

    def integrate_boundary_layer(self, n):
        """ Integrates the boundary-layer and calculates the scale factor from displacement thickness
        arg: int: n: Iteration number from the iterative solver.
        returns: ndarray: y: Wall normal coordinates.
        returns: ndarray: u: Streamwise velocity component profile.
        returns: ndarray: T: Temperature profile.
        returns: float: scale: Scale factor of the boundary-layer.
        """
        sumd, record_z = 0, 0
        z = np.zeros(n+1)
        d_eta = self.eta[1]*0.5
        # self.soln[1,:] is the u velocity, should be 1 in free stream
        for i in range(1, n+1):
            z[i] = z[i-1] + d_eta*(self.soln[3, i] + self.soln[3, i-1])
            dm = self.soln[3, i-1] - self.soln[1, i-1]
            dd = self.soln[3, i] - self.soln[1, i]
            sumd += d_eta*(dd+dm)
            if(self.soln[1, i] > 0.999 and record_z < 1.0):
                # print "recording at iteration: ", i
                dlta = z[i]
                record_z = 2.0
            scale = sumd
        print "delta is :", dlta
        print "conversion factor is: ", scale
        print "scaled delta is: ", dlta/scale
        # Rescale with displacement thickness and convert to FLOWER variable normalisation
        y, u, T = z/scale, self.soln[1, :], self.soln[3, :]
        return y, u, T, scale


class Initialise_Katzer(object):

    def __init__(self, npoints, lengths, directions, beta, n_coeffs, block, coordinate_strings):
        """ Generates the initialiastion equations for the boundary-layer profile.
        arg: list: npoints: Numerical values of the number of points in each direction.
        arg: list: lengths: Numerical values of the problem dimensions.
        arg: list: directions: Integer values of the problem directions.
        arg: float: beta: Stretching factor for stretched grids.
        arg: int: n_coeffs: Desired number of coefficients for the polynomial fit.
        arg: object: block: OpenSBLI SimulationBlock.
        """
        self.block = block
        self.directions = directions
        self.idxs = block.grid_indexes
        self.npoints = npoints
        self.lengths = lengths
        self.betas = beta
        self.n_coeffs = n_coeffs
        self.Tinf = 1.0
        self.Tw = 1.67619431
        self.coordinate_strings = coordinate_strings
        self.initial = self.generate_initial_condition()
        return

    def generate_initial_condition(self):
        # Load u, T profile
        y, u, T, rho, n = self.load_compbl()
        n_coeffs = self.n_coeffs
        # y, u, T, rho, n = load_similarity(115.0, 255, 5.0)
        tolerance = 1e-9
        npoints, lengths, betas = self.npoints, self.lengths, self.betas
        if len(self.directions) == 1:  # 2D Katzer and 3D spanwise periodic Katzer
            direction = 1
            # Uniform [0, 1] values to perform grid stretching
            eta = np.linspace(0, 1, npoints[direction])
            coordinates = self.generate_coordinates(npoints, lengths, betas, eta, y)
            x, y_stretch = coordinates[0], coordinates[1]
            # Interpolate u, T, rho onto the grid
            u_new = self.interpolate_onto_grid(y, y_stretch, u, 0.51425, 0)
            T_new = self.interpolate_onto_grid(y, y_stretch, T, 0, 0)
            rho_new = 1.0/T_new
            rhou_new = rho_new*u_new
            # Solve continuity equation to obtain rhov
            rhov_new = self.solve_continuity(x, y_stretch, u_new, rho_new)
            edge = self.find_edge_of_bl(u_new, tolerance)
            # Obtain polynomial fit coefficients
            rhou_coeffs = self.fit_polynomial(y_stretch, rhou_new, edge, n_coeffs)
            rhov_coeffs = self.fit_polynomial(y_stretch, rhov_new, edge, n_coeffs)
            T_coeffs = self.fit_polynomial(y_stretch, T_new, edge, n_coeffs)
            self.generate_one_wall_equations([rhou_new, rhov_new, T_new], [rhou_coeffs, rhov_coeffs, T_coeffs], direction, edge)
        elif len(self.directions) == 2:  # 3D with one side wall in x2
            directions = [1, 2]
            z = np.linspace(0, lengths[2], npoints[2])
            edges, coeffs, profiles, normal_coeffs, normal_profiles = [], [], [], [], []
            for dire in directions:
                # Uniform [0, 1] values to perform grid stretching
                eta = np.linspace(0, 1, npoints[dire])
                coordinates = self.generate_coordinates(npoints, lengths, betas, eta, y, z)
                x, y_stretch, z_stretch = coordinates[0], coordinates[1], coordinates[2]
                # Interpolate u, and T onto the grid
                u_new = self.interpolate_onto_grid(y, coordinates[dire], u, 0.51425, 0)
                T_new = self.interpolate_onto_grid(y, coordinates[dire], T, 0, 0)
                # Temperature scaling function in region [0, 1]
                g = self.temperature_scaling(T_new)
                rho_new = 1.0/T_new
                rhou_new = rho_new*u_new
                # Solve continuity equation to obtain rho*wall_normal_velocity
                rho_vel_normal = self.solve_continuity(x, coordinates[dire], u_new, rho_new)
                profiles.append([rhou_new, T_new])
                normal_profiles.append(rho_vel_normal)
                edge = self.find_edge_of_bl(u_new, tolerance)
                edges.append(edge)
                # Obtain polynomial fit coefficients
                rhou_coeffs = self.fit_polynomial(coordinates[dire], rhou_new, edge, n_coeffs)
                g_coeffs = self.fit_polynomial(coordinates[dire], g, edge, n_coeffs)
                rho_vel_normal_coeffs = self.fit_polynomial(coordinates[dire], rho_vel_normal, edge, n_coeffs)
                coeffs.append([rhou_coeffs, g_coeffs])
                normal_coeffs.append(rho_vel_normal_coeffs)
            self.generate_two_wall_equations(profiles, coeffs, directions, edges, normal_profiles, normal_coeffs)
        else:
            raise NotImplementedError("Boundary layer initialisation is only implemented for 1 or two walls.")
        # Output the initialisation
        initial = GridBasedInitialisation()
        initial.add_equations(self.coordinate_strings + self.eqns)
        return initial

    def temperature_scaling(self, temp_profile):
        """ Computes the temperature profile between [0,1]
        args: ndarray: temp_profile: Temperature profile values between wall temperature and freestream.
        returns: ndarray: g: Temperature profile ranging from 0 to 1. """
        g = (temp_profile - self.Tw)/(self.Tinf - self.Tw)
        return g

    def form_equation(self, variable, name, coefficients, direction, edge):
        """ Creates the piecewise equations for the cases of 2D and 3D span-periodic boundary-layer profiles.
        arg: ndarray: variable: Array of values for a given flow variable, used to obtain the free-stream value.
        arg: string: name: Name of the variable.
        arg: ndarray: coefficients: Coefficients for the polynomial fit.
        arg: int: direction: Spatial direction to apply the equation to.
        arg: int: edge: Grid index for the edge of the boundary-layer.
        returns: Eq: eqn: OpenSBLI equation to add to the initialisation kernel.
        """
        powers = [i for i in range(np.size(coefficients))][::-1]
        eqn = sum([coeff*GridVariable('x1')**power for (coeff, power) in zip(coefficients, powers)])
        idx = self.idxs[direction]
        eqn = Eq(GridVariable('%s' % name), Piecewise((eqn, idx < edge), (variable[edge], True)))
        return eqn

    def form_mixed_equation(self, profiles, names, coefficients, directions, edges, normal_profiles, normal_coeffs):
        """ Generates the equations for the 3D SBLI sidewall case.
        arg: list: profiles: Arrays of values for the rhou and [0,1] temperature profiles.
        arg: list: names: Variable names as strings.
        arg: list: coefficients: Coefficients for the polynomial fit for rhou and temperature profiles.
        arg: list: directions: Directions that contain a wall.
        arg: list: edges: Indices for the boundary layer edges in each direction.
        arg: list: normal_profiles: Arrays of values for the rhov and rhow profiles.
        arg: list: normal_coeffs: Coefficients for the polynomial fit for rhov and rhow.
        returns: list: piecewise_eqns: Piecewise initialisation equations to be added to the initialisation class."""
        idx0, idx1 = self.idxs[directions[0]], self.idxs[directions[1]]
        # Assuming we have the same number of poly coefficients in each direction, change later if required
        powers = [i for i in range(np.size(coefficients[0][0]))][::-1]
        # Loop over rhou, and T profiles
        piecewise_eqns = []
        # Create rhou profiles
        eqn1 = sum([coeff*GridVariable('x%d' % directions[0])**power for (coeff, power) in zip(coefficients[0][0], powers)])
        eqn2 = sum([coeff*GridVariable('x%d' % directions[1])**power for (coeff, power) in zip(coefficients[1][0], powers)])  # profiles[0][max(edges)]
        u_var1, u_var2 = GridVariable('%s_1' % names[0]), GridVariable('%s_2' % names[0])
        freestream_value = np.max([profiles[0][0][edges[0]], profiles[1][0][edges[1]]])
        freestream_value = 1.0
        piecewise_eqns.append(Eq(u_var1, Piecewise((eqn1, idx0 < edges[0]), (freestream_value, True))))
        piecewise_eqns.append(Eq(u_var2, Piecewise((eqn2, idx1 < edges[1]), (freestream_value, True))))
        piecewise_eqns.append(Eq(GridVariable('%s' % names[0]), u_var1*u_var2))
        # Create T profiles, g = (T-Tw)/(Tinf - Tw) ---> T = g*(Tinf-Tw) + Tw
        Tw = ConstantObject('Twall')
        eqn1 = sum([coeff*GridVariable('x%d' % directions[0])**power for (coeff, power) in zip(coefficients[0][1], powers)])
        eqn2 = sum([coeff*GridVariable('x%d' % directions[1])**power for (coeff, power) in zip(coefficients[1][1], powers)])  # profiles[0][max(edges)]
        T_var1, T_var2 = GridVariable('%s_1' % names[1]), GridVariable('%s_2' % names[1])
        freestream_value = np.max([profiles[0][1][edges[0]], profiles[1][1][edges[1]]])
        freestream_value = 1.0
        piecewise_eqns.append(Eq(T_var1, Piecewise((eqn1, idx0 < edges[0]), (freestream_value, True))))
        piecewise_eqns.append(Eq(T_var2, Piecewise((eqn2, idx1 < edges[1]), (freestream_value, True))))
        piecewise_eqns.append(Eq(GridVariable('%s' % names[1]), T_var1*T_var2*(freestream_value - Tw) + Tw))
        # Create normal velocity component profiles, rhov:
        Lx2 = ConstantObject('Lx2')
        eqn1 = sum([coeff*GridVariable('x%d' % directions[0])**power for (coeff, power) in zip(normal_coeffs[0], powers)])
        temp1 = GridVariable('%s' % names[2])
        rhov_inf = normal_profiles[0][edges[0]]*u_var2  # Multiplying by rhou from the other direction
        piecewise_eqns.append(Eq(temp1, Piecewise((eqn1*u_var2, idx0 < edges[0]), (rhov_inf, True))))
        # rhow:
        eqn2 = sum([coeff*GridVariable('x%d' % directions[1])**power for (coeff, power) in zip(normal_coeffs[1], powers)])
        temp2 = GridVariable('%s' % names[3])
        # rhow should reduce to zero at the symmetry plane
        rhow_inf = normal_profiles[1][edges[1]]*(1 - GridVariable('x2')/Lx2)*u_var1  #Multiplying by rhou from the other direction
        piecewise_eqns.append(Eq(temp2, Piecewise((eqn2*u_var1, idx1 < edges[1]), (rhow_inf, True))))
        return piecewise_eqns

    def generate_one_wall_equations(self, data, coeffs, direction, edge):
        """ Generates the equations for 2D SBLI and 3D span-periodic cases.
        arg: list: data: Profile arrays for rhou0, rhou1 and temperature.
        arg: list: coeffs: Coefficients for the polynomial fits.
        arg: list: direction: Direction normal to the wall.
        arg: int: edge: Grid index for the edge of the boundary-layer."""
        self.eqns = [Eq(GridVariable('x%d' % direction), DataObject('x%d' % direction))]
        rhou0_eqn = self.form_equation(data[0], 'rhou0', coeffs[0], direction, edge)
        rhou1_eqn = self.form_equation(data[1], 'rhou1', coeffs[1], direction, edge)
        T_eqn = self.form_equation(data[2], 'T', coeffs[2], direction, edge)
        # Set conservative values
        rho, rhou0, rhou1, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('T')
        rho_eqn = Eq(rho, 1.0/T)
        rho_store = Eq(DataObject('rho'), rho)
        rhou0_store = Eq(DataObject('rhou0'), rhou0)
        rhou1_store = Eq(DataObject('rhou1'), rhou1)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = Eq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2)/rho)
        self.eqns += [rhou0_eqn, rhou1_eqn, T_eqn, rho_eqn, rho_store, rhou0_store, rhou1_store, rhoE_store]
        if self.block.ndim == 3:  # Periodic case, rhow = 0
            self.eqns += [Eq(DataObject('rhou2'), 0.0)]
        return

    def generate_two_wall_equations(self, profiles, coeffs, directions, edges, normal_profile, normal_coeffs):
        """ Generates the equations for 2D SBLI and 3D span-periodic cases.
        arg: list: profiles: Profile arrays for rhou0 and temperature.
        arg: list: coeffs: Coefficients for the polynomial fits.
        arg: list: directions: Directions normal to the wall.
        arg: list: edges: Grid indexes for the edge of the boundary-layer.
        arg: list: normal_profile: Profile for the wall normal velocity components.
        arg: list: normal_coeffs: Coefficients for the wall normal polynomial fit."""
        # Create GridVariables of the boundary-layer independent variables x1, x2
        self.eqns = [Eq(GridVariable('x%d' % direction), DataObject('x%d' % direction)) for direction in directions]
        names = ['rhou0', 'T', 'rhou1', 'rhou2']
        # Create the piecewise equations formed from the boundary layers
        bl_equations = self.form_mixed_equation(profiles, names, coeffs, directions, edges, normal_profile, normal_coeffs)
        # rhou0, T, rhou1, rhou2 = bl_equations[0], bl_equations[1], bl_equations[2], bl_equations[3]
        # Set conservative values
        rho, rhou0, rhou1, rhou2, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('rhou2'), GridVariable('T')
        rho_eqn = Eq(rho, 1.0/T)
        rho_store = Eq(DataObject('rho'), rho)
        rhou0_store = Eq(DataObject('rhou0'), rhou0)
        rhou1_store = Eq(DataObject('rhou1'), rhou1)
        rhou2_store = Eq(DataObject('rhou2'), rhou2)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        rhoE_store = Eq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2 + rhou2**2)/rho)
        self.eqns += bl_equations + [rho_eqn, rho_store, rhou0_store, rhou1_store, rhou2_store, rhoE_store]
        return

    def function_to_fit(self, x, a, b, c, d):
        """
        General form of an equation to fit to the data, to obtain the coefficients.
        """
        # return a*np.exp(-c*(x-b))+d
        return a*x**3 + b*x**2 + c*x + d

    def fit_polynomial(self, coords, variable, bl_edge, n_coeffs):
        """ Fits a polynomial to the input data, coefficients are returned.
        arg: ndarray: coords: Independent variable of the input data.
        arg: ndarray: variable: Dependent variable of the input data.
        arg: int: bl_edge: Array index at the edge of the boundary-layer.
        arg: int: n_coeffs: Desired number of coefficients for the polynomial.
        returns: ndarray: coeffs: Coefficients of the polynomial fit.
        """
        coords = coords[0:bl_edge]
        variable = variable[0:bl_edge]
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                coeffs = poly.polyfit(coords, variable, n_coeffs)
            except np.RankWarning:
                print "Poorly conditioned fit"
        ffit = poly.polyval(coords, coeffs)
        plt.plot(coords, ffit, label='fit')
        plt.plot(coords, variable, label='original_data')
        plt.legend(loc="best")
        # plt.savefig("domain_fit.png", dpi=200)
        # plt.show()
        # Reverse coefficients so they are in descending order
        return coeffs[::-1]

    def fit_function(self, coords, variable, bl_edge, func):
        coords = coords[0:bl_edge]
        variable = variable[0:bl_edge]
        popt, pcov = curve_fit(func, coords, variable)
        plt.plot(coords, func(coords, *popt), 'r-', label='fit')
        plt.plot(coords, variable, 'k-', label='original_data')
        plt.legend()
        # plt.show()
        # print popt, pcov
        return

    def solve_continuity(self, x, y, u, rho):
        """
        Solves the continuity equation to obtain the wall normal velocity profile.
        arg: ndarray: x: Independent coordinate values.
        arg: ndarray: y: Dependent coordinate values.
        arg: ndarray: u: Streamwise velocity component values.
        arg: ndarray: rho: Density values.
        returns: ndarray: rhov: Array of values for the wall normal velocity components.
        """
        # Grid offset delta to form derivative approximation
        n = np.size(y)
        ya2 = y[:]
        delta, scale, re = 0.1, 2.316725, 950.0
        rex0 = 0.5*(re/scale)**2
        x0 = 0.5*re/scale**2
        drudx, rhov = np.zeros_like(y), np.zeros_like(y)
        # Local Reynolds number scaling to obtain a v profile
        sqrex = np.sqrt(rex0 + re*x[0])
        delsx = np.sqrt(2.0)*scale*(x0+x[0])/sqrex
        ya2 = delsx*ya2
        d2y_u = spline(ya2, u, n, 0.51425, 0)
        d2y_rho = spline(ya2, rho, n, 0, 0)
        for j in range(0, n):
            ya2[j] = delsx*ya2[j]
            dstarp = delsx*np.sqrt((x0+x[0]+delta)/(x0+x[0]))
            dstarm = delsx*np.sqrt((x0+x[0]-delta)/(x0+x[0]))
            yp, ym = ya2[j]/dstarp, ya2[j]/dstarm
            uxp = splint(ya2, u, d2y_u, n, yp)
            uxm = splint(ya2, u, d2y_u, n, ym)
            rhoxp = splint(ya2, rho, d2y_rho, n, yp)
            rhoxm = splint(ya2, rho, d2y_rho, n, ym)
            drudx[j] = (rhoxp*uxp-rhoxm*uxm)/(2.0*delta)
            rhov[j] = rhov[j-1]-0.5*(ya2[j]-ya2[j-1])*(drudx[j]+drudx[j-1])
        return rhov

    def load_similarity(self, length, npoints, beta):
        bl = Boundary_layer_profile(2.0, 0.72, 1.4, -1, 950, length, npoints, beta)
        y, u, T, rho, n = bl.y, bl.u, bl.T, 1.0/bl.T, np.size(bl.y)
        return y, u, T, rho, n

    def load_compbl(self):
        data = np.loadtxt("./comp_bl.dat")
        y, u, T = data[:, 0], data[:, 1], data[:, 2]
        rho, n = 1.0/T, np.size(y)
        return y, u, T, rho, n

    def generate_coordinates(self, npoints, lengths, betas, eta, y, z=None):
        """ Creates the coordinate arrays.
        arg: list: npoints: Number of points in each direction.
        arg: list: lengths: Domain length in each direction.
        arg: list: betas: Stretch factors in each direction.
        arg: ndarray: eta: Uniform values between [0,1] used for stretching.
        arg: ndarray: y: Coordinate values in the y direction.
        arg: ndarray: z: Coordinate values in the z direction.
        returns: list: coordinates: List of arrays containing the coordinate values."""
        x = np.linspace(0, lengths[0], npoints[0])
        if self.block.ndim == 2:
            by = betas[0]
            y_stretch = y[-1]*np.sinh(by*eta)/np.sinh(by)
            coordinates = [x, y_stretch]
        elif self.block.ndim == 3:
            by, bz = betas[0], betas[1]
            y_stretch = y[-1]*np.sinh(by*eta)/np.sinh(by)
            z_stretch = z[-1]*np.sinh(bz*eta)/np.sinh(bz)
            coordinates = [x, y_stretch, z_stretch]
        return coordinates

    def generate_uniform_coordinates(self, npoints, lengths, betas, y, z=None):
        coordinates = []
        for i in range(len(npoints)):
            coordinates.append(np.linspace(0, lengths[i], npoints[i]))
        return coordinates

    def interpolate_onto_grid(self, y_in, y_out, var_in, y0, yn):
        # Create interpolating second derivative spline
        # n = size of the original data
        n = np.size(y_in)
        n_out = np.size(y_out)
        d2y = spline(y_in, var_in, n, y0, yn)
        # Array for variable interpolated onto the grid
        var_out = np.zeros_like(y_out)
        for i in range(n_out):
            var_out[i] = splint(y_in, var_in, d2y, n, y_out[i])
        return var_out

    def find_edge_of_bl(self, variable, tolerance):
        """ Finds the edge of the boundary layer and returns the index of that grid point.
        arg: ndarray: variable: Array of values for a given flow variable.
        arg: float: tolerance: Stopping tolerance for the difference between two successive grid points.
        returns: int: index: Index of the boundary-layer edge."""
        index = 1
        while np.abs(variable[index]-variable[index-1]) > tolerance:
            index += 1
        return index

    def freestream_value(self, variable, index):
        """ Returns the value of a flow variable for a given grid index.
        arg: ndarray: variable: Array of values for a given flow variable.
        arg: index: Array index to use.
        returns: float: variable[index]: Value of the flow variable at that index."""
        return variable[index]
