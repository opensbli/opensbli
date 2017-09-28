import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py
import glob
import sys
import os.path
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from sympy import *
import matplotlib.cm as cm
import matplotlib.transforms as transforms


plt.style.use('classic')

# Matplotlib settings for publication-ready figures
try:
    f = open(os.path.expanduser('./rcparams.py'), 'r')
    exec(f.read())
except:
    pass

def read_file(fname):
    # Read in the simulation output
    dump = glob.glob("./" + fname)
    if not dump or len(dump) > 1:
        print "Error: No dump file found, or more than one dump file found."
        sys.exit(1)
    f = h5py.File(dump[-1], 'r')
    group = f["opensbliblock00"]
    return f, group

def extract_data(group, lhalo, rhalo, k):

    #linear dimensions of the dataset
    np = group["rho_B0"].shape
    rho = group["rho_B0"].value
    rhou = group["rhou0_B0"].value
    rhoE = group["rhoE_B0"].value

    rho = rho[lhalo:np[0]-rhalo]
    grid_points = [np[0] - 2*k]
    rhou = rhou[lhalo:np[0]-rhalo]
    rhoE = rhoE[lhalo:np[0]-rhalo]

    u = rhou/rho
    p = (0.4)*(rhoE - 0.5*(u**2)*rho)
    a = numpy.sqrt(1.4*p/rho)
    M = u/a
    T = 1.4*4*p/rho
    return rho, u, rhoE, p, M, T

def line_graphs(x, variable, name):
    plt.plot(x, variable)
    plt.xlabel(r'$x_0$', fontsize=20)
    plt.ylabel(r'$%s$' % name, fontsize=20)
    plt.savefig("%s.pdf" % name, bbox_inches='tight')
    plt.clf()
    return

def plot(fname, n_levels):
    f, group1 = read_file(fname)
    rho, u, rhoE, P, M, T = extract_data(group1, 5, 5, 3)
    x = numpy.linspace(0, 10, 1600)
    coordinates = [x]
    variables = [rho, u, rhoE, P, M, T]
    names = ["\\rho", "u", "\\rho E", "P", "M", "T"]

    # Line plots1
    variables = [rho[:], u[:], P[:]/P[0]]
    names = ["\\rho", "u", "P"]
    for var, name in zip(variables, names):
        line_graphs(x, var, name)
    plt.clf()
    f.close()

fname = "opensbli.h5"
plot(fname, 25)
