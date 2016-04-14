""" Reads in the 2D MMS solution data and computes the error between the numerical and analytical/manufactured solution. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import string
from scipy.interpolate import griddata
import glob
import sys

# Matplotlib settings for publication-ready figures
try:
    f = open('/home/ctj1r15/opensbli-paper/images/rcparams.py', 'r')
    exec(f.read())
except:
    pass

fig, subax = plt.subplots(2, 3, sharex='col', sharey='row', figsize=(9,5))
subax = subax.flat
    
def plot_phi(simulation_index, phi, phi_analytical):
    """ Plot the field 'phi' in a multi-figure plot, to show the field converging to the analytical solution. """
    subax[simulation_index].imshow(phi, extent=(0, 2*pi, 0, 2*pi), interpolation='nearest', aspect='auto')
    
    if simulation_index != 0 and simulation_index != 1 and simulation_index != 2:
        subax[simulation_index].set_xlabel(r'$x$')
    subax[simulation_index].set_ylabel(r'$y$')
    subax[simulation_index].text(0.5, 1.05, "("+string.ascii_lowercase[simulation_index]+")", transform=subax[simulation_index].transAxes, size=10, weight='bold')
    if simulation_index == 4: # The last simulation.
        # Add in the analytical solution as well as a special case.
        im = subax[5].imshow(phi_analytical, extent=(0, 2*pi, 0, 2*pi), aspect='auto')
        subax[5].set_xlabel(r'$x$')
        subax[5].set_ylabel(r'$y$')
        
        # Insert colourbar
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax, label=r"$\phi$")
        
        subax[5].text(0.5, 1.05, "("+string.ascii_lowercase[5]+")", transform=subax[5].transAxes, size=10, weight='bold')
        
        plt.savefig("phi.pdf", bbox_inches='tight')
        plt.clf()
    
    
    
def compute_error(degree, simulation_index, number_of_points):
    # Number of grid points. This is assumed to be the same in the x and y directions.
    nx = number_of_points
    ny = number_of_points  
    
    # Number of halo nodes at each end
    halo = degree/2

    # Read in the simulation output
    path = "./mms_%d_%d/mms_%d_%d_opsc_code/" % (degree, simulation_index, degree, simulation_index)
    dump = glob.glob(path + "/mms_*.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["mms_%d_%d_block" % (degree, simulation_index)]
    
    # Get the numerical solution field
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at the left (and bottom) end of the domain. Include one strip of halo points at the right (and top) of the domain to enforce periodicity in the solution field plot.
    phi = phi[halo:nx+halo+1, halo:ny+halo+1]
    print phi.shape
    
    # Grid spacing. Note: The length of the domain is divided by nx (or ny) and not nx-1 (or ny-1) because of the periodicity. In total we have nx+1 points, but we only solve nx points; the (nx+1)-th point is set to the same value as the 0-th point to give a full period, to save computational effort.
    dx = (2.0*pi)/(nx)
    dy = (2.0*pi)/(ny)
    
    # Coordinate arrays. 
    x = numpy.zeros((nx+1)*(ny+1)).reshape((nx+1, ny+1))
    y = numpy.zeros((nx+1)*(ny+1)).reshape((nx+1, ny+1))

    phi_analytical = numpy.zeros((nx+1)*(ny+1)).reshape((nx+1, ny+1))

    # Compute the error by first interpolating the numerical and analytical results onto a much finer grid of points and computing the L2 norm of the absolute difference.
    grid_points = []
    grid_numerical = []
    grid_analytical = []
    target_grid_x, target_grid_y = numpy.mgrid[0:2*pi:32j, 0:2*pi:32j]
    for i in range(0, nx+1):
        for j in range(0, ny+1):
            # Work out the x and y coordinates. Note the swapping of the 'j' and 'i' here.
            x[i,j] = j*dx
            y[i,j] = i*dy
            grid_points.append([x[i,j], y[i,j]])
            grid_numerical.append(phi[i,j])
            phi_analytical[i,j] = sin(x[i,j])*cos(y[i,j])
            grid_analytical.append(phi_analytical[i,j])

    grid_points = numpy.array(grid_points)
    grid_numerical = numpy.array(grid_numerical)
    grid_analytical = numpy.array(grid_analytical)
    interpolated_numerical = griddata(grid_points, grid_numerical, (target_grid_x, target_grid_y), method='nearest')
    interpolated_analytical = griddata(grid_points, grid_analytical, (target_grid_x, target_grid_y), method='nearest')

    # Only plot phi for the 6th order simulations.
    if degree == 12:
        plot_phi(simulation_index, phi, phi_analytical)
    
    return numpy.linalg.norm(abs(interpolated_numerical - interpolated_analytical), ord=2)

def plot():
    # Plot the error against the grid spacing dx.
    Lx = 2*pi
    
    degrees = range(2, 13, 2)
    errors = []
    dxs = []
    for d in range(len(degrees)):
        dx = []
        error = []
        for simulation_index in range(0, 5):
            number_of_points = 4*(2**simulation_index)
            dx.append(Lx/number_of_points)            
            error.append(compute_error(degrees[d], simulation_index, number_of_points))
        errors.append(error)
        dxs.append(dx)
    
    print "Errors in the L2 norm: ", errors
    plt.clf()
    
    colours_expected = ["-r", "-g", "-b", "-y", "-c", "-k"]
    colours = ["o--r", "o--g", "o--b", "o--y", "o--c", "o--k"]
    for d in range(0, len(degrees)):
        # Plot the errors.
        plt.loglog(dxs[d], errors[d], colours[d], label=r"Order = %d" % degrees[d])
        
        # Plot the expected convergence line for comparison.
        expected_convergence = (numpy.array([0.8*max(errors[d])*(1.0/2**degrees[d])**i for i in range(len(errors[d]))]))
        plt.loglog(dxs[d], expected_convergence, colours_expected[d], label=r"")
        
    plt.xlabel(r"Grid spacing $\Delta x$ (m)")
    plt.ylabel(r"Solution error in the L2 norm")
    plt.legend(loc='best', numpoints=1)
    plt.savefig("mms_convergence_analysis.pdf", bbox_inches='tight')
    
    
if(__name__ == "__main__"):
    plot()
