""" Reads in the 2D MMS solution data, written by the OPS dat writer, and plots the fields. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def plot(path):
    # Number of grid points
    nx = 10
    ny = 10
    
    nu = 1.0
    
    # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["mms_block"]
    
    u = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    u = u[halo:nx+halo, halo:ny+halo]
    print u.shape
    print u

    
    # Grid spacing
    L = 1.0
    dx = (2.0*pi*L)/(nx)
    dy = (2.0*pi*L)/(ny)
    
    # Coordinate arrays
    x = numpy.zeros(nx*ny).reshape((nx, ny))
    y = numpy.zeros(nx*ny).reshape((nx, ny))

    # Compute the error
    t = 10.0
    for i in range(0, nx):
        for j in range(0, ny):
            x[i,j] = i*dx
            y[i,j] = j*dy

    plt.imshow(u)
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$\phi$")
    plt.legend()
    plt.savefig("phi.pdf")
    


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
