""" Reads in the 1D wave equation solution data, written by the OPS dat writer,
and plots the scalar field 'phi'. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def plot(path):
    # Number of grid points
    nx = 8
    
     # Number of halo nodes at each end
    halo = 2
    
    # Time-step size
    t = 0.1

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["wave_block"]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo]
    print phi
    # Grid spacing
    dx = 1/(nx);
    
    # Coordinate array
    x = numpy.zeros(nx)

    # Compute the error
    for i in range(0, nx):
        x[i] = i*dx
        # Analytical solution
        #u_analytical[i,j] = cos(x[i,j])*sin(y[i,j])*exp(-2*nu*t)
        #v_analytical[i,j] = -sin(x[i,j])*cos(y[i,j])*exp(-2*nu*t)
        #u_error[i,j] = u[i,j] - u_analytical[i,j]
        #v_error[i,j] = v[i,j] - v_analytical[i,j]

    plt.imshow(phi)
    plt.show()


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
