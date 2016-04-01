""" Reads in the 2D MMS solution data and computes the error between the numerical and analytical/manufactured solution. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import string

fig, subax = plt.subplots(2, 3, sharex='col', sharey='row', figsize=(9,5))
subax = subax.flat

def plot_phi(simulation_index, phi, phi_analytical):
    """ Plot the field 'phi' in a multi-figure plot, to show the field converging to the analytical solution. """
    subax[simulation_index].imshow(phi, extent=(0, 10, 0, 10), interpolation='nearest', aspect='auto')
    
    if simulation_index != 0 and simulation_index != 1 and simulation_index != 2:
        subax[simulation_index].set_xlabel(r'$x$')
    subax[simulation_index].set_ylabel(r'$y$')
    subax[simulation_index].text(0.5, 1.05, "("+string.ascii_lowercase[simulation_index]+")", transform=subax[simulation_index].transAxes, size=10, weight='bold')
    if simulation_index == 4: # The last simulation.
        # Add in the analytical solution as well as a special case.
        im = subax[5].imshow(phi_analytical, extent=(0, 10, 0, 10), aspect='auto')
        subax[5].set_xlabel(r'$x$')
        subax[5].set_ylabel(r'$y$')
        
        # Insert colourbar
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax, label=r"$\phi$")
        
        subax[5].text(0.5, 1.05, "("+string.ascii_lowercase[5]+")", transform=subax[5].transAxes, size=10, weight='bold')
        
        plt.savefig("phi.pdf", bbox_inches='tight')
        plt.clf()
    
    
    
def error(simulation_index, number_of_points):
    # Number of grid points. This is assumed to be the same in the x and y directions.
    nx = number_of_points
    ny = number_of_points  
    
    # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    path = "./mms_%d/mms_%d_opsc_code/" % (simulation_index, simulation_index)
    f = h5py.File(path + "/state.h5", 'r')
    group = f["mms_%d_block" % simulation_index]
    
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
            phi_analytical[i,j] = sin(x[i,j])*cos(y[i,j])
            phi_error[i,j] = abs(phi[i,j] - phi_analytical[i,j])

    plot_phi(simulation_index, phi, phi_analytical)
    
    return numpy.linalg.norm(phi_error, ord=2)

def plot():
    # Plot the error against the grid spacing dx.
    Lx = 2*pi
    dx = []
    errors = []
    for simulation_index in range(0, 5):
        number_of_points = 4*(2**simulation_index)
        dx.append(Lx/number_of_points)
        errors.append(error(simulation_index, number_of_points))
    print "Errors in the L2 norm: ", errors
    plt.loglog(dx, errors, 'o-k', label=r"$\phi$")
    
    # Plot the third-order convergence line for comparison.
    third_order = (numpy.array([1e-1*(1.0/8.0)**i for i in range(len(errors))]))
    plt.loglog(dx, third_order, '--k', label=r"Third-order convergence")
    
    plt.xlabel(r"Grid spacing $\Delta x$ (m)")
    plt.ylabel(r"Error in the L2 norm")
    plt.legend(loc='best', numpoints=1)
    plt.savefig("mms_convergence_analysis.pdf", bbox_inches='tight')
    
    
if(__name__ == "__main__"):
    plot()
