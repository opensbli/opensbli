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
    
    # Simulation parameters
    nu = 1.0/1600.0 # 1/Re
    t = 1000*0.0005
    gamma = 1.4

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["auto_block_OPSC[0]"]
    
    u_im = group["rhou0[0]"].value
    v_im = group["rhou1[0]"].value
    r = group["rho[0]"].value
    E = group["rhoE[0]"].value

    # Obtain the velocity field by dividing the momentum through by the density.
    u_im /= r
    v_im /= r
    
    # Ignore the 2 halo nodes at either end of the domain
    u = u_im[halo:nx+halo, halo:ny+halo]
    v = v_im[halo:nx+halo, halo:ny+halo]

    # Grid spacing
    dx = 2.0*pi/(nx);
    dy = 2.0*pi/(ny);
    
    # Coordinate arrays
    x = numpy.zeros(nx*ny).reshape((nx, ny))
    y = numpy.zeros(nx*ny).reshape((nx, ny))
    
    # Analytical solution arrays
    u_analytical = numpy.zeros(nx*ny).reshape((nx, ny))
    v_analytical = numpy.zeros(nx*ny).reshape((nx, ny))
    
    # Error arrays
    u_error = numpy.zeros(nx*ny).reshape((nx, ny))
    v_error = numpy.zeros(nx*ny).reshape((nx, ny))

    # Compute the error
    for i in range(0, nx):
        for j in range(0, ny):
            x[i,j] = i*dx
            y[i,j] = j*dy
            # Analytical solution
            u_analytical[i,j] = cos(x[i,j])*sin(y[i,j])*exp(-2*nu*t)
            v_analytical[i,j] = -sin(x[i,j])*cos(y[i,j])*exp(-2*nu*t)
            u_error[i,j] = u[i,j] - u_analytical[i,j]
            v_error[i,j] = v[i,j] - v_analytical[i,j]

    print u_analytical
    plt.imshow(u)
    plt.show()


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="tgv_2d_error")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
