""" Reads in the 2D Taylor-Green vortex solution data, written by the OPS dat writer,
and computes the error between the numerical and analytical solution. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def plot(path):
    # Number of grid points in the x and y directions
    nx = 32
    ny = 32
    halo = 2 # Number of halo nodes at each end

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["auto_block_OPSC[0]"]
    
    c = group["c[0]"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    c = c[halo:nx+halo, halo:ny+halo]

    # Grid spacing
    dx = 2.0*pi/(nx);
    dy = 2.0*pi/(ny);
    
    # Coordinate arrays
    x = numpy.zeros(nx*ny).reshape((nx, ny))
    y = numpy.zeros(nx*ny).reshape((nx, ny))
    print c
    plt.imshow(c)
    plt.show()


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="tgv_2d_error")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
