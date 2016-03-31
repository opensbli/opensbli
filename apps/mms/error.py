""" Reads in the 2D MMS solution data and computes the error between the numerical and analytical/manufactured solution. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def error(path, i, number_of_points):
    # Number of grid points. This is assumed to be the same in the x and y directions.
    nx = number_of_points
    ny = number_of_points
    
    # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["mms_%d_block" % i]
    
    # Get the numerical solution field
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo, halo:ny+halo]
    print phi.shape
    
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
    print phi
    plt.imshow(phi)
    plt.show()
    return numpy.linalg.norm(phi_error, ord=2)

def plot(path):
    
    errors = []
    
    for i in range(0, 5):
        number_of_points = 10*(2**i)
        errors.append(error(path + "/mms_%d/mms_%d_opsc_code/" % (i, i), i, number_of_points))
    
    print errors

if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
