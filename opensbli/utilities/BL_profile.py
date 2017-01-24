import numpy as np
from sympy import *
from scipy import interpolate

class GridGen(object):
	def __init__(self):
		self.generate_stretched_x1(255, 115.0, 5.0)
		return

	def generate_stretched_x1(self, np_x1, L_x1, stretch_factor):
		""" Generates a grid stretched in the wall normal (x1) direction. """
		b_x1 = stretch_factor
		# Eta is [0,1] uniformly distributed over the number of points in y direction
		eta = np.linspace(0,1,np_x1)
		x1 = L_x1*np.sinh(b_x1*eta)/np.sinh(b_x1)
		np.savetxt('py_grid.dat', x1, fmt='%.7e', delimiter=' ', newline='\n', header=str(np_x1))
		return


class Boundary_layer_profile(object):
	def __init__(self):
		self.generate_boundary_layer_profile(None, 2.0, 0.72, 1.4, -1.0, 950.0)
		return

	def compbl(self, v):
		"""Specifies the ODE for the boundary-layer."""
		# Sutherland law
		suth = self.suth
		c = np.sqrt(v[3])*(1.0+suth)/(v[3]+suth)
		dcdg = 1.0/(2.0*np.sqrt(v[3])) - np.sqrt(v[3])/(v[3]+suth)
		dcdg *= (1.0+suth) / (v[3]+suth)
		cp = dcdg*v[4]

		dv = np.zeros(5)
		dv[0] = v[1]
		dv[1] = v[2]
		dv[2] = -v[2]*(cp+v[0])/c
		dv[3] = v[4]
		dv[4] = -v[4]*(cp+self.pr*v[0])/c - self.pr*(self.gama-1)*self.xmach**2 * v[2]**2
		return dv

	def RK4(self, y, dydx, n, h):
		""" Runge-Kutta routine. """
		nmax = 5
		y_out = np.zeros(n)
		yt, dyt, dym = np.zeros(nmax), np.zeros(nmax), np.zeros(nmax)
		hh, h6 = h*0.5, h/6.0
		yt = y + hh*dydx
		dyt = self.compbl(yt)
		yt = y + hh*dyt
		dym = self.compbl(yt)
		yt = y + h*dym
		dym = dyt + dym
		dyt = self.compbl(yt)
		y_out = y + h6*(dydx + dyt + 2.0*dym)
		return y_out

	def RKDUMB(self, vstart, n_vars, x1, x2, nstep):
		nmax = 5
		v, dv = np.zeros(nmax), np.zeros(nmax)
		v = vstart[:]
		self.y[:,0] = v[:]
		self.xx[0] = x1
		x = x1
		h = (x2-x1)/float(nstep)
		for k in range(0, nstep):
			dv = self.compbl(v)
			v = self.RK4(v, dv, n_vars, h)
			x += h
			self.xx[k+1] = x
			self.y[:,k+1] = v[:]
		return

	def generate_boundary_layer_profile(self, y_points, xmach, pr, gama, Tw, Re):
		""" Generates a boundary layer initial profile. Solves the mean flow
		in a compressible boundary layer. (Equations 7.32) in White (1974).
		arg txtfile: Text file containing stretched y-coordinates.
		arg float: xmach: Mach number.
		arg float: pr: Prandtl number.
		arg float: gama: Ratio of specific heats.
		arg float: Tw: Wall temperature Tw/Tinf (< 0 for adiabatic)
		arg float: Re: Reynolds number."""
		# Declaring arrays
		nmax, maxit, jmax = 5, 10, 1001
		self.xx, self.y = np.zeros(jmax), np.zeros((nmax,jmax))
		v, dv, f, f1, f2 = np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)

		self.pr, self.gama, self.xmach, self.Re, self.Tw = pr, gama, xmach, Re, Tw
		self.nvisc = 1 # 1 Sutherlands law, 2 Power law, 3 Chapman-Rubesin approximation.
		Tinf = 288.0
		etamax = 10.0
		nstep = jmax-1 
		errtol = 1e-10
		
		self.suth = 110.40/Tinf # Sutherland constant 
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
		for k in range(0,5):
			vstart[2] = v[0] # Initial value
			if Tw < 0:
				vstart[3] = v[1]
				twall = vstart[3]
			else:
				vstart[4] = 1.141575 # v[1]
			# Calling RKDUMB here
			start = time.time()
			self.RKDUMB(vstart, 5, 0.0, etamax, nstep)
			end = time.time()
			# print "Time for RKDUMB was : ", (end-start)
			err = np.abs(self.y[1,nstep]-1.0) + np.abs(self.y[3,nstep]-1.0)
			# print "The error is: ", err
			# print "y_out last value is ", self.y[:, nstep]
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
			self.RKDUMB(vstart, 5, 0.0, etamax, nstep)
			# print "y_out last value is ", self.y[:, nstep]
			f1[0] = self.y[1, nstep] - 1.0
			f1[1] = self.y[3, nstep] - 1.0

			# Increment v[1]
			vstart[2] = v[0]
			if (Tw < 0):
				vstart[3] = v[1] + dv[1]
				twall = vstart[3]
			else:
				vstart[4] = v[1] + dv[1]
			self.RKDUMB(vstart, 5, 0.0, etamax, nstep)
			# print "y_out last value is ", self.y[:, nstep]
			f2[0] = self.y[1, nstep] - 1.0
			f2[1] = self.y[3, nstep] - 1.0
			# print "dv values: ", dv
			# Solve linear system
			al11 = (f1[0] - f[0])/dv[0]
			al21 = (f1[1] - f[1])/dv[0]
			al12 = (f2[0] - f[0])/dv[1]
			al22 = (f2[1] - f[1])/dv[1]
			det = float(al11*al22 - al21*al12)
			# New dv for improved solution
			dv[0] = (-al22*f[0] + al12*f[1])/det
			dv[1] = (al21*f[0] - al11*f[1])/det
			v[0] = v[0] + dv[0]
			v[1] = v[1] + dv[1]
			# Write out improved estimate
			print 'it   ', k, dv[0], dv[1], v[0], v[1]
		print "final vaues ", v[0], v[1], err
		print "calling output with nstep: ", nstep
		self.output(nstep)
		return

	def output(self, n):
		sumd, record_z = 0, 0
		yy, z = np.zeros(n+1), np.zeros(n+1)
		u2, t2 = np.zeros(n+1), np.zeros(n+1)
		## self.y[1,:] is u velocity, should be 1 in free stream
		## if condition ends when at 0.999 of freestream velocity
		for i in range(1,n+1):
			z[i] = z[i-1] + 0.5*(self.xx[i] - self.xx[i-1])*(self.y[3,i] + self.y[3,i-1])
			dm = self.y[3,i-1] - self.y[1,i-1]
			dd = self.y[3,i] - self.y[1,i]
			sumd += 0.5*(self.xx[i] - self.xx[i-1])*(dd+dm)
			if(self.y[1,i] > 0.999 and record_z < 1.0):
				print "recording at iteration: ", i
				dlta = z[i]
				record_z = 2.0
			scale = sumd
		print "delta is :", dlta
		print "conversion factor is: ", scale
		print "scaled delta is: ", dlta/scale
		# Rescale with displacement thickness and convert to FLOWER variable normalisation
		y, u, T = z/scale, self.y[1,:], self.y[3,:]
		# Call SPLINE routine with 1st derivatives zero at the start/end points
		u2 = interpolate.CubicSpline(y, u, bc_type='clamped')
		t2 = interpolate.CubicSpline(y, T, bc_type='clamped')
		# Read in the stretched grid coordinates
		yy = np.loadtxt('./py_grid.dat', skiprows=1)
		uu, tt = u2(yy), t2(yy)
		# Find the index where the u0 velocity becomes 1.0 
		# and set the u0, T entries to 1.0 from this point to match SBLI
		idx = np.argmin(np.abs(uu - 1.0))
		uu[idx+2:], tt[idx+2:] = 1.0, 1.0
 		return yy, uu, tt

	def generate_2D_grid(nx, ny, y_points):
		""" Takes the stretched y coordinates, and generates the initial condition in
		rho, rhou0, rhou1, rhoE, plus the x, and y grid coordinates. """
		return
