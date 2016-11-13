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

""" Reads in the 1D viscous Burgers solution data and plots the scalar field 'phi'. """

import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys

def plot(path):
    # Number of grid points
    nx = 200
    
     # Number of halo nodes at each end
    halo = 2

    # Read in the simulation output
    dump = glob.glob(path + "/viscous_burgers_*.h5")
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["viscous_burgers_block"]
    
    phi = group["phi"].value
    
    # Ignore the 2 halo nodes at either end of the domain
    phi = phi[halo:nx+halo]
    print phi
    # Grid spacing
    dx = 1.0/(nx);
    
    # Coordinate array
    x = numpy.zeros(nx)
    phi_analytical = numpy.zeros(nx)
    phi_error = numpy.zeros(nx)

    # Compute the error
    for i in range(0, nx):
        x[i] = i*dx

    plt.plot(x, phi, label=r"Solution field $\phi$")
    plt.xlabel(r"$x$ (m)")
    plt.ylabel(r"$\phi$")
    plt.legend()
    plt.savefig("phi.pdf", bbox_inches='tight')
    


if(__name__ == "__main__"):
    plot(path="./viscous_burgers_opsc_code")
