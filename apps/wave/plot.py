#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

""" Reads in the 1D wave equation solution data, written by the OPS dat writer,
and plots the scalar field 'phi'. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys
import os.path

# Matplotlib settings for publication-ready figures
try:
    f = open(os.path.expanduser('~/opensbli-paper/images/rcparams.py'), 'r')
    exec(f.read())
except:
    pass
    
def plot(path):
    # Number of grid points
    nx = 1000
    
     # Number of halo nodes at each end
    halo = 4

    # Read in the simulation output
    dump = glob.glob(path + "/wave_*.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["wave_block"]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo]

    # Grid spacing
    dx = 1.0/(nx);
    
    # Coordinate array
    x = numpy.zeros(nx)
    phi_initial = numpy.zeros(nx)
    phi_analytical = numpy.zeros(nx)
    phi_error = numpy.zeros(nx)

    # Compute the error
    for i in range(0, nx):
        x[i] = i*dx
        # Initial condition
        phi_initial[i] = sin(2*pi*x[i])
        # Analytical solution
        phi_analytical[i] = sin(2*pi*(x[i]+0.5)) # Phi should be a sin wave shifted to the right by x = 0.5 (since the wave speed is 0.5 m/s and we've simulated until T = 1.0 s).
        phi_error[i] = abs(phi_analytical[i] - phi[i])
    
    plt.clf()

    plt.plot(x, phi_error, "-k", label=r"Absolute error in $\phi(x,\ t=1)$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"Error")
    plt.legend()
    plt.savefig("phi_error.pdf", bbox_inches='tight')

    plt.clf()
    plt.plot(x, phi_initial, "--k", label=r"$\phi(x,\ t=0)$")
    plt.plot(x, phi, "-k", label=r"$\phi(x,\ t=1)$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"Wave amplitude $\phi$")
    plt.legend()
    plt.savefig("phi.pdf", bbox_inches='tight')
    
    plt.clf()    


if(__name__ == "__main__"):
    plot("./wave_opsc_code")
