""" Reads in the 2D MMS solution data and computes the error between the numerical and analytical/manufactured solution. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def error(i, number_of_points):
    # Number of grid points. This is assumed to be the same in the x and y directions.
    nx = number_of_points
    ny = number_of_points  
    
    # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    path = "./mms_%d/mms_%d_opsc_code/" % (i, i)
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
            # Work out the x and y coordinates. Note the swapping of the 'j' and 'i' here.
            x[i,j] = j*dx
            y[i,j] = i*dy
            phi_analytical[i,j] = sin(x[i,j])
            phi_error[i,j] = abs(phi[i,j] - phi_analytical[i,j])
    return numpy.linalg.norm(phi_error, ord=2)

def plot():
    # Plot the error against the grid spacing dx.
    Lx = 2*pi
    dx = []
    errors = []
    for i in range(0, 5):
        number_of_points = 4*(2**i)
        dx.append(Lx/number_of_points)
        errors.append(error(i, number_of_points))
    print "Errors in the L2 norm: ", errors
    plt.loglog(dx, errors, '-k', label=r"$\phi$")
    
    # Plot the third-order convergence line for comparison.
    third_order = (numpy.array([1e-1*(1.0/8.0)**i for i in range(len(errors))]))
    plt.loglog(dx, third_order, '--k', label=r"Third-order convergence")
    
    plt.xlabel(r"Grid spacing $\Delta x$ (m)")
    plt.ylabel(r"Error in the L2 norm")
    plt.legend(loc='best')
    plt.show()
    
    
if(__name__ == "__main__"):
    plot()
