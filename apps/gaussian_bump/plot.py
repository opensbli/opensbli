#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

""" Plots the solution field phi from the 2D advection-diffusion simulation. This should be a Gaussian bump advected in the positive x-direction. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys

def plot(path):
    # Number of grid points in the x and y directions
    nx = 10
    ny = 10
    halo = 2 # Number of halo nodes at each end

    # Read in the simulation output
    dump = glob.glob(path + "/gaussian_bump_*.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["gaussian_bump_block"]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo, halo:ny+halo]

    # Grid spacing
    L = 10.0
    dx = L/nx
    dy = L/ny
    
    # Coordinate arrays
    x = numpy.zeros(nx*ny).reshape((nx, ny))
    y = numpy.zeros(nx*ny).reshape((nx, ny))
    print phi
    plt.imshow(phi)
    plt.show()


if(__name__ == "__main__"):
    plot("./gaussian_bump_opsc_code")
