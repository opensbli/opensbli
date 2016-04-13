""" Reads in the 1D wave equation solution data, written by the OPS dat writer,
and plots the scalar field 'phi'. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def plot(path, simulation_index):
    # Number of grid points
    nx = 10*(2**simulation_index)
    
     # Number of halo nodes at each end
    halo = 4

    # Read in the simulation output
    f = h5py.File(path + "/state.h5", 'r')
    group = f["wave_%d_block" % simulation_index]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo]

    # Grid spacing
    dx = 1.0/(nx);
    
    # Coordinate array
    x = numpy.zeros(nx)
    phi_initial = numpy.zeros(nx)
    phi_analytical = numpy.zeros(nx)
    phi_error = numpy.zeros(nx)

    # Compute the error
    from scipy.interpolate import griddata
    grid_x = numpy.mgrid[0:1:10000j]
    for i in range(0, nx):
        x[i] = i*dx
        # Initial condition
        phi_initial[i] = sin(2*pi*x[i])
        # Analytical solution
        phi_analytical[i] = sin(2*pi*(x[i]+0.5)) # Phi should be a sin wave shifted to the right by x = 0.5 (since the wave speed is 0.5 m/s and we've simulated until T = 1.0 s).
        phi_error[i] = abs(phi_analytical[i] - phi[i])
    
    data1 = griddata(x, phi, grid_x, method='nearest')
    data2 = griddata(x, phi_analytical, grid_x, method='nearest')
    #print data
    print numpy.linalg.norm(abs(data1-data2), ord=2)
    
    plt.clf()

    plt.plot(x, phi_error, "-k", label=r"Absolute error in $\phi(x,\ t=1)$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"Error")
    plt.legend()
    plt.savefig("phi_error.pdf", bbox_inches='tight')

    plt.clf()
    plt.plot(x, phi_initial, "--k", label=r"$\phi(x,\ t=0)$")
    plt.plot(x, phi, "-k", label=r"$\phi(x,\ t=1)$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"Wave amplitude $\phi$")
    plt.legend()
    plt.savefig("phi.pdf", bbox_inches='tight')
    
    plt.clf()    


if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the directory containing the output files.", action="store", type=str)
    args = parser.parse_args()
    
    for simulation_index in range(0, 5):
        plot("wave_%d/wave_%d_opsc_code" % (simulation_index, simulation_index), simulation_index)
