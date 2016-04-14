""" Plots the solution field phi from the 2D advection-diffusion simulation. This should be a Gaussian bump advected in the positive x-direction. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys

def plot(path):
    # Number of grid points in the x and y directions
    nx = 10
    ny = 10
    halo = 2 # Number of halo nodes at each end

    # Read in the simulation output
    dump = glob.glob(path + "/gaussian_bump_*.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["gaussian_bump_block"]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo, halo:ny+halo]

    # Grid spacing
    L = 10.0
    dx = L/nx
    dy = L/ny
    
    # Coordinate arrays
    x = numpy.zeros(nx*ny).reshape((nx, ny))
    y = numpy.zeros(nx*ny).reshape((nx, ny))
    print phi
    plt.imshow(phi)
    plt.show()


if(__name__ == "__main__"):
    plot("./gaussian_bump_opsc_code")
