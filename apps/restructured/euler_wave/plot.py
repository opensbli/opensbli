import numpy
import matplotlib.pyplot as plt
import h5py
import glob
import sys
import os.path
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib.transforms as transforms


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
    f = h5py.File(fname, 'r')
    group = f["opensbliblock00"]
    return f, group

def read_dataset(group, dataset):
    d_m = group["%s" % (dataset)].attrs['d_m']
    size = group["%s" % (dataset)].shape
    read_start = [abs(d) for d in d_m]
    read_end = [s-abs(d) for d, s in zip(d_m, size)]
    if len(read_end) == 2:
        read_data = group["%s" % (dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1]]
    elif len(read_end) == 3:
        read_data = group["%s" % (dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
    else:
        raise NotImplementedError("")
    return read_data

def extract_data(group):
    # for key in group.keys():
    #     data_without_halos = read_dataset(group, key)
    #     output_opensbli.create_dataset("%s" % (key), data=data_without_halos)
    rho = read_dataset(group, "rho_B0")
    rhou = read_dataset(group, "rhou0_B0")
    rhov = read_dataset(group, "rhou1_B0")
    rhoE = read_dataset(group, "rhoE_B0")
    x = read_dataset(group, "x0_B0")
    y = read_dataset(group, "x1_B0")
    u = rhou/rho
    v = rhov/rho
    return x, y, rho, u, v, rhoE

def plot(fname, n_levels):
    f, group1 = read_file(fname)
    x, y, rho, u, v, rhoE = extract_data(group1)

    coordinates = [x, y]
    variables = [rho, u, v, rhoE]
    names = ["\\rho", "u", "v", "\\rho E"]

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

    # Constant values for exact solution
    v_const = -0.5
    u_const = 1.0
    t = 2.5
    exact = 1.0 + 0.2*numpy.sin(numpy.pi*(x+y - t*(u_const+v_const)))
    npoints = numpy.shape(exact)
    print numpy.shape(exact)

    rho_error = numpy.abs(exact - rho)
    L1 = numpy.sum(rho_error)/(npoints[0]*npoints[1])
    Linf = numpy.max(rho_error)
    errors = [L1, Linf]

    text_file = open("errors.txt", "w")
    text_file.write("L1, Linf\n")
    text_file.write("%e, %e"%(L1, Linf))
    text_file.close()
    print "=================================="
    print "L^1 error: %e " % L1
    print "L_inf error: %e " % Linf
    f.close()

fname = "opensbli_output.h5"
plot(fname, 5)
