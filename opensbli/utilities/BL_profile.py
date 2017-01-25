import numpy as np
from sympy import *
from scipy.interpolate import CubicSpline
from scipy.integrate import odeint

class GridGen(object):
    def __init__(self):
        return

    def generate_stretched_x1(self, np_x1, L_x1, stretch_factor):
        """ Generates a grid stretched in the wall normal (x1) direction. """
        b_x1 = stretch_factor
        # Eta is [0,1] uniformly distributed over the number of points in y direction
        eta = np.linspace(0,1,np_x1)
        x1 = L_x1*np.sinh(b_x1*eta)/np.sinh(b_x1)
        return x1

class Boundary_layer_profile(object):
    def __init__(self):
        gr = GridGen()
        y_points = gr.generate_stretched_x1(255, 115.0, 5.0)
        # y_points = np.linspace(0, 115.0, 255)
        u, T = self.generate_boundary_layer_profile(y_points, 2.0, 0.72, 1.4, -1.0, 950.0)
        for i in range(np.shape(u)[0]):
            print y_points[i], u[i], T[i]
        return

    def compbl(self, v, p=None):
        suth = self.suth
        c = np.sqrt(v[3])*(1.0+suth)/(v[3]+suth)
        dcdg = 1.0/(2.0*np.sqrt(v[3])) - np.sqrt(v[3])/(v[3]+suth)
        dcdg *= (1.0+suth) / (v[3]+suth)
        cp = dcdg*v[4]
        dv = [v[1], v[2], -v[2]*(cp+v[0])/c, v[4], \
             -v[4]*(cp+self.pr*v[0])/c - self.pr*(self.gama-1)*self.xmach**2 * v[2]**2]
        return dv

    def generate_boundary_layer_profile(self, y_points, xmach, pr, gama, Tw, Re):
        """ Generates a boundary layer initial profile. Solves the mean flow
        in a compressible boundary layer. (Equations 7.32) in White (1974).
        arg txtfile: Text file containing stretched y-coordinates.
        arg float: xmach: Mach number.
        arg float: pr: Prandtl number.
        arg float: gama: Ratio of specific heats.
        arg float: Tw: Wall temperature Tw/Tinf (< 0 for adiabatic)
        arg float: Re: Reynolds number."""
        nmax, maxit, jmax = 5, 10, 1001
        v, dv, f, f1, f2 = np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)
        self.pr, self.gama, self.xmach, self.Re, self.Tw = pr, gama, xmach, Re, Tw
        self.nvisc = 1 # 1 Sutherlands law, 2 Power law, 3 Chapman-Rubesin approximation.
        Tinf = 288.0
        etamax = 10.0
        nstep = jmax-1 
        errtol = 1e-10
        self.suth = 110.40/Tinf # Sutherland constant 
        self.y = np.zeros((nmax,jmax))
        self.xx = np.linspace(0, etamax, jmax)
        vstart = np.zeros(5)
        if Tw > 0:
            vstart[3] = Tw
        elif Tw < 0: # Adiabatic wall, initial guess
            v[1] = 1.0 + 0.5*0.84*(gama-1)*xmach**2
            v[0] = 0.47*v[1]**0.21
        else: # Fixed wall temperature, initial guess. Tested for 1<Tw<4 and M<8 (Pr=0.72, gamma=1.4)
            v[1] = 4.5 #!0.062*xmach**2-0.1*(Tw-1.0)*(1.0+xmach)/(0.2+xmach)
            v[0] = 0.494891 #!0.45-0.01*xmach+(Tw-1.)*0.06
            twall = Tw
        # Initial increments
        dv[0] = v[0]*0.01
        if Tw < 0:
            dv[1] = v[1]*0.01
        else:
            dv[1] = 0.1
        # Main loop
        for k in range(nmax):
            vstart[2] = v[0] # Initial value
            if Tw < 0:
                vstart[3] = v[1]
                twall = vstart[3]
            else:
                vstart[4] = 1.141575 # v[1]
            # ODE solver parameters
            abserr = 1.0e-10
            relerr = 1.0e-8
            # Call the ODE solver.
            self.y = odeint(self.compbl, vstart, self.xx,
                     atol=abserr, rtol=relerr).T
            err = np.abs(self.y[1,nstep]-1.0) + np.abs(self.y[3,nstep]-1.0)
            if err < errtol:
                break
            f[0] = self.y[1, nstep] - 1.0
            f[1] = self.y[3, nstep] - 1.0
            # Increment v[0]
            vstart[2] = v[0] + dv[0]
            if (Tw < 0):
                vstart[3] = v[1]
                twall = vstart[3]
            else:
                vstart[4] = v[1]
            self.y = odeint(self.compbl, vstart, self.xx,
              atol=abserr, rtol=relerr).T
            f1[0] = self.y[1, nstep] - 1.0
            f1[1] = self.y[3, nstep] - 1.0
            # Increment v[1]
            vstart[2] = v[0]
            if (Tw < 0):
                vstart[3] = v[1] + dv[1]
                twall = vstart[3]
            else:
                vstart[4] = v[1] + dv[1]
            self.y = odeint(self.compbl, vstart, self.xx,
              atol=abserr, rtol=relerr).T
            f2[0] = self.y[1, nstep] - 1.0
            f2[1] = self.y[3, nstep] - 1.0
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
        u, T = self.output(y_points, nstep)
        return u, T

    def output(self, y_points, n):
        sumd, record_z = 0, 0
        z = np.zeros(n+1)
        # self.y[1,:] is the u velocity, should be 1 in free stream
        for i in range(1,n+1):
            z[i] = z[i-1] + 0.5*(self.xx[i] - self.xx[i-1])*(self.y[3,i] + self.y[3,i-1])
            dm = self.y[3,i-1] - self.y[1,i-1]
            dd = self.y[3,i] - self.y[1,i]
            sumd += 0.5*(self.xx[i] - self.xx[i-1])*(dd+dm)
            if(self.y[1,i] > 0.999 and record_z < 1.0):
                # print "recording at iteration: ", i
                dlta = z[i]
                record_z = 2.0
            scale = sumd
        print "delta is :", dlta
        print "conversion factor is: ", scale
        print "scaled delta is: ", dlta/scale
        # Rescale with displacement thickness and convert to FLOWER variable normalisation
        y, u, T = z/scale, self.y[1,:], self.y[3,:]
        # Call the SPLINE routine with 1st derivatives zero at the start/end points
        u2 = CubicSpline(y, u, bc_type='clamped')
        T2 = CubicSpline(y, T, bc_type='clamped')
        # Perform the interpolation to the stretched y coordinates
        u, T = u2(y_points), T2(y_points)
        # Find the index where the u0 velocity becomes 1.0 
        # and set the u0, T entries to 1.0 from this point to match SBLI
        idx = np.argmin(np.abs(u - 1.0))
        u[idx+2:], T[idx+2:] = 1.0, 1.0
        return u, T

    def generate_2D_grid(nx, ny, y_points):
        """ Takes the stretched y coordinates, and generates the initial condition in
        rho, rhou0, rhou1, rhoE, plus the x, and y grid coordinates. """
        return