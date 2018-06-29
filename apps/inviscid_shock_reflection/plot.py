import numpy
import matplotlib.pyplot as plt
import h5py
import glob
import sys
import os.path
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from sympy import *
import matplotlib.cm as cm


plt.style.use('classic')


# Matplotlib settings for publication-ready figures
try:
    f = open(os.path.expanduser('./rcparams.py'), 'r')
    exec(f.read())
except:
    pass


def contour_local(fig, levels0, label, x, y, variable):
    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    ax1.set_xlabel(r"$x_0$", fontsize=20)
    ax1.set_ylabel(r"$x_1$", fontsize=20)
    CS = ax1.contourf(x, y, variable, levels=levels0, cmap=cm.jet)
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    ticks_at = numpy.linspace(levels0[0], levels0[-1], 10)
    cbar = plt.colorbar(CS, cax=cax1, ticks=ticks_at, format='%.3f')
    cbar.ax.set_ylabel(r"$%s$" % label, fontsize=20)
    return


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
    # linear dimensions of the dataset
    np = group["rho_B0"].shape
    rho = group["rho_B0"].value
    rhou = group["rhou0_B0"].value
    rhov = group["rhou1_B0"].value
    rhoE = group["rhoE_B0"].value
    x = group["x0_B0"].value
    y = group["x1_B0"].value

    rho = rho[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    grid_points = [np[0] - 2*k, np[1] - 2*k]
    rhou = rhou[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    rhov = rhov[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    x = x[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    y = y[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    rhoE = rhoE[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]

    u = rhou/rho
    v = rhov/rho
    p = (0.4)*(rhoE - 0.5*(u**2+v**2)*rho)
    a = numpy.sqrt(1.4*p/rho)
    M = u/a
    T = 1.4*4*p/rho
    return x, y, rho, u, v, rhoE, p, M, T


def line_graphs(x, variable, name):
    if name == "u":
        plt.axhline(y=0.0, linestyle='--', color='k')

    plt.plot(x[1, :], variable)
    plt.xlabel(r'$x_0$', fontsize=20)
    plt.ylabel(r'$%s$ at wall' % name, fontsize=20)
    plt.savefig("wall_%s.pdf" % name, bbox_inches='tight')
    plt.clf()
    return


def plot(fname, n_levels):
    f, group1 = read_file(fname)
    x, y, rho, u, v, rhoE, P, M, T = extract_data(group1, 5, 5, 3)
    coordinates = [x, y]
    variables = [rho, u, v, rhoE, P, M, T]
    names = ["\\rho", "u", "v", "\\rho E", "P", "M", "T"]

    # Contour plots
    for var, name in zip(variables, names):
        min_val = numpy.min(var)
        max_val = numpy.max(var)
        levels = numpy.linspace(min_val, max_val, n_levels)
        print "%s" % name
        print levels
        fig = plt.figure()
        contour_local(fig, levels, "%s" % name, x, y, var)
        plt.savefig("fig_%s.pdf" % name, bbox_inches='tight')
        plt.clf()

    # Line plots1
    variables = [rho[0, :], u[1, :], P[0, :]/P[0, 0]]
    names = ["\\rho", "u", "P"]
    for var, name in zip(variables, names):
        line_graphs(x, var, name)
    # Inlet temperature profile
    plt.semilogy(T[:, 0], y[:, 0])
    plt.ylabel('x1')
    plt.xlabel('T at inlet')
    plt.savefig("temperature.pdf", bbox_inches='tight')
    plt.clf()
    # V velocity at top boundary
    plt.plot(x[1, :], v[-2, :])
    plt.xlabel('x0')
    plt.ylabel('V velocity at top boundary')
    plt.savefig("V_top.pdf", bbox_inches='tight')
    plt.clf()
    f.close()


fname = "opensbli_output.h5"
plot(fname, 25)
