""" Reads in the 2D MMS solution data and computes the error between the numerical and analytical/manufactured solution. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def error(path):
    # Number of grid points
    nx = 10
    ny = 10
    
    # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["mms_block"]
    
    # Get the numerical solution field
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo, halo:ny+halo]
    print phi.shape
    print phi

    # Grid spacing
    dx = (2.0*pi)/(nx)
    dy = (2.0*pi)/(ny)
    
    # Coordinate arrays
    x = numpy.zeros(nx*ny).reshape((nx, ny))
    y = numpy.zeros(nx*ny).reshape((nx, ny))

    phi_analytical = numpy.zeros(nx*ny).reshape((nx, ny))
    phi_error = numpy.zeros(nx*ny).reshape((nx, ny))
    
    # Compute the error
    for i in range(0, nx):
        for j in range(0, ny):
            x[i,j] = i*dx
            y[i,j] = j*dy
            phi_analytical[i,j] = sin(x[i,j])
            phi_error[i,j] = abs(phi[i,j] - phi_analytical[i,j])
    
    return numpy.linalg.norm(phi_error, order=2)


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    print "The error is %f" % error(args.path)
