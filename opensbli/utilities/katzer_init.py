from opensbli.utilities.numerical_functions import spline, splint
from sympy import Eq, Piecewise
from scipy.integrate import odeint
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from opensbli.initialisation import GridBasedInitialisation
from opensbli.core.kernel import Kernel
from opensbli.core.opensbliobjects import DataObject, ConstantObject
from opensbli.core.grid import GridVariable

plt.style.use('classic')


class Boundary_layer_profile(object):

    def __init__(self, xmach, Pr, gama, Tw, Re, length, npoints, beta):
        self.y, self.u, self.T, self.scale = self.generate_boundary_layer_profile(xmach, Pr, gama, Tw, Re)
        self.Re = Re
        self.interpolate(length, npoints, beta)
        self.n = np.size(self.y)
        return

    def interpolate(self, Ly, npoints, beta):
        eta = np.linspace(0, 1, npoints)
        # y_full = np.linspace(0, Ly, npoints)
        y_full = Ly*np.sinh(beta*eta)/np.sinh(beta)
        y, u, T = self.y[:], self.u[:], self.T[:]
        n = self.n
        d2y_u = spline(y, u, n, 0, 0)
        d2y_T = spline(y, T, n, 0, 0)
        u_full, T_full = np.zeros_like(y_full), np.zeros_like(y_full)
        n_full = np.size(u_full)
        for i in range(n_full):
            # Stop spline interpolation if at the edge of the boundary-layer
            if (y_full[i] < y[-1]):
                u_full[i] = splint(y, u, d2y_u, n, y_full[i])
                T_full[i] = splint(y, T, d2y_T, n, y_full[i])
            else:
                u_full[i] = u_full[i-1]
                T_full[i] = T_full[i-1]
        self.y, self.u, self.T = y_full, u_full, T_full
        return

    def compbl(self, v, p=None):
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

    def __init__(self, npoints, lengths, directions, beta, n_coeffs, block):
        self.block = block
        self.directions = directions
        self.idxs = block.grid_indexes
        self.npoints = npoints
        self.lengths = lengths
        self.beta = beta
        self.n_coeffs = n_coeffs
        return

    def main(self):
        # Load u, T profile
        y, u, T, rho, n = self.load_compbl()
        # y, u, T, rho, n = load_similarity(115.0, 255, 5.0)
        x, y_stretch = self.generate_coordinates(self.npoints, self.lengths, y, self.beta)
        # Interpolate u, T, rho onto the grid
        u_new = self.interpolate_onto_grid(y, y_stretch, u, 0.51425, 0)
        T_new = self.interpolate_onto_grid(y, y_stretch, T, 0, 0)
        rho_new = 1.0/T_new
        rhou_new = rho_new*u_new
        # Solve continuity equation to obtain rhov
        rhov_new = self.solve_continuity(x, y_stretch, u_new, rho_new)
        edge = self.find_edge_of_bl(u_new)
        # Obtain polynomial fit coefficients
        n_coeffs = self.n_coeffs
        if len(self.directions) == 1:
            direction = 1
            rhou_coeffs = self.generate_poly_coefficients(y_stretch, rhou_new, edge, n_coeffs)
            rhov_coeffs = self.generate_poly_coefficients(y_stretch, rhov_new, edge, n_coeffs)
            T_coeffs = self.generate_poly_coefficients(y_stretch, T_new, edge, n_coeffs)
            self.generate_2D_equations([rhou_new, rhov_new, T_new], [rhou_coeffs, rhov_coeffs, T_coeffs], direction, edge)
        elif len(self.directions) == 2:
            pass
        else:
            raise NotImplementedError("Boundary layer initialisation is only implemented for 1 or two walls.")

        return

    def form_equation(self, variable, name, coefficients, direction, edge):
        powers = [i for i in range(np.size(coefficients))][::-1]
        eqn = sum([coeff*GridVariable('x1')**power for (coeff, power) in zip(coefficients, powers)])
        idx = self.idxs[direction]
        eqn = Eq(GridVariable('%s' % name), Piecewise((eqn, idx < edge), (variable[edge], True)))
        return eqn

    def generate_2D_equations(self, data, coeffs, direction, edge):
        self.eqns = [Eq(GridVariable('x%d' % direction), DataObject('x%d' % direction))]
        eqn1 = self.form_equation(data[0], 'rhou0', coeffs[0], direction, edge)
        eqn2 = self.form_equation(data[1], 'rhou1', coeffs[1], direction, edge)
        eqn3 = self.form_equation(data[2], 'T', coeffs[2], direction, edge)
        # Set conservative values
        rho, rhou0, rhou1, T = GridVariable('rho'), GridVariable('rhou0'), GridVariable('rhou1'), GridVariable('T')
        eqn4 = Eq(rho, 1.0/T)
        eqn5 = Eq(DataObject('rho'), rho)
        eqn6 = Eq(DataObject('rhou0'), rhou0)
        eqn7 = Eq(DataObject('rhou1'), rhou1)
        gama, Minf = ConstantObject("gama"), ConstantObject("Minf")
        eqn8 = Eq(DataObject('rhoE'), rho*T/(gama*(gama-1)*Minf**2) + 0.5*(rhou0**2 + rhou1**2)/rho)
        self.eqns += [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8]
        return

    def create_initial_condition(self):
        initial = GridBasedInitialisation()

        initial.add_equations([])
        return initial

    def function_to_fit(self, x, a, b, c, d):
        # return a*np.exp(-c*(x-b))+d
        return a*x**3 + b*x**2 + c*x + d

    def fit_polynomial(self, coords, variable, bl_edge, n_coeffs):
        coords = coords[0:bl_edge]
        variable = variable[0:bl_edge]
        coeffs = poly.polyfit(coords, variable, n_coeffs)
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
        print popt, pcov
        return

    def solve_continuity(self, x, y_stretch, u, rho):
        # Grid offset delta to form derivative approximation
        n = np.size(y_stretch)
        ya2 = y_stretch[:]
        delta, scale, re = 0.1, 2.316725, 950.0
        rex0 = 0.5*(re/scale)**2
        x0 = 0.5*re/scale**2
        drudx, rhov = np.zeros_like(y_stretch), np.zeros_like(y_stretch)
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

    def generate_coordinates(self, npoints, lengths, y, beta):
        # Uniform in x for now
        x = np.linspace(0, lengths[0], npoints[0])
        eta = np.linspace(0, 1, npoints[1])
        y_stretch = y[-1]*np.sinh(beta*eta)/np.sinh(beta)
        return x, y_stretch

    def interpolate_onto_grid(self, y_in, y_out, var_in, y0, yn):
        # Create interpolating second derivative spline
        # n = size of the original data
        n = np.size(y_in)
        d2y = spline(y_in, var_in, n, y0, yn)
        # Array for variable interpolated onto the grid
        var_out = np.zeros_like(y_out)
        for i in range(n):
            var_out[i] = splint(y_in, var_in, d2y, n, y_out[i])
        return var_out

    def generate_poly_coefficients(self, coordinates, variable, end_of_bl, n_coeffs):
        coeffs = self.fit_polynomial(coordinates, variable, end_of_bl, n_coeffs)
        return coeffs

    def find_edge_of_bl(self, variable):
        i = 1
        while np.abs(variable[i]-variable[i-1]) > 1e-16:
            i += 1
        print "edge is :", i
        return i

    def freestream_value(self, variable, index):
        return variable[index]
