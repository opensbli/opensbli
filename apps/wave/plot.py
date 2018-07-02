""" Reads in the 1D wave equation solution data, written by the OPS dat writer,
and plots the scalar field 'phi'. """

import numpy
from math import pi, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys
import os.path


def plot(path):
    # Number of grid points
    nx = 200

    # Number of halo nodes at each end
    halo = 5

    # Read in the simulation output
    dump = glob.glob(path + "/opensbli_output.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["opensbliblock00"]

    phi = group["phi_B0"].value
    x = group["x0_B0"].value

    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo]
    x = x[halo:nx+halo]

    # Grid spacing
    dx = 1.0/(nx)

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
        phi_analytical[i] = sin(2*pi*(x[i]+0.5))  # Phi should be a sin wave shifted to the right by x = 0.5 (since the wave speed is 0.5 m/s and we've simulated until T = 1.0 s).
        phi_error[i] = abs(phi_analytical[i] - phi[i])

    plt.clf()

    plt.plot(x, phi_error, "-k", label=r"Absolute error in $\phi(x,\ t=1)$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"Error")
    plt.legend()
    plt.savefig(directory + "phi_error.pdf", bbox_inches='tight')

    plt.clf()
    plt.plot(x, phi_initial, "--k", label=r"$\phi(x,\ t=0)$")
    plt.plot(x, phi, "-k", label=r"$\phi(x,\ t=1)$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"Wave amplitude $\phi$")
    plt.legend()
    plt.savefig(directory + "phi.pdf", bbox_inches='tight')
    plt.clf()


directory = './simulation_plots/'

if not os.path.exists(directory):
    os.makedirs(directory)

plot('.')
