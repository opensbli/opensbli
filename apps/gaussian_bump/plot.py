""" Plots the solution field phi from the 2D advection-diffusion simulation. This should be a Gaussian bump advected in the positive x-direction. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def plot(path):
    # Number of grid points in the x and y directions
    nx = 10
    ny = 10
    halo = 2 # Number of halo nodes at each end

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
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
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
