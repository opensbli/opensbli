""" Reads in the 1D viscous Burgers solution data and plots the scalar field 'phi'. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys

def plot(path):
    # Number of grid points
    nx = 200
    
     # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    dump = glob.glob(path + "/viscous_burgers_*.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["viscous_burgers_block"]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo]
    print phi
    # Grid spacing
    dx = 1.0/(nx);
    
    # Coordinate array
    x = numpy.zeros(nx)
    phi_analytical = numpy.zeros(nx)
    phi_error = numpy.zeros(nx)

    # Compute the error
    for i in range(0, nx):
        x[i] = i*dx

    plt.plot(x, phi, label=r"Solution field $\phi$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$\phi$")
    plt.legend()
    plt.savefig("phi.pdf", bbox_inches='tight')
    


if(__name__ == "__main__"):
    plot(path="./viscous_burgers_opsc_code")
