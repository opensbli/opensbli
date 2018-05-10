import numpy
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm


plt.style.use('classic')


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


def plot(files, n_levels, min_val, max_val):
    levels = numpy.linspace(min_val, max_val, n_levels)

    for num, file in files:
        print "Processing image: %d" % num
        f, group = read_file(file)
        np = group["rho_B0"].shape
        rho = group["rho_B0"].value
        rho = rho[5:-5, 5:-5]
        x, y = group["x0_B0"].value, group["x1_B0"].value
        x, y = x[5:-5, 5:-5], y[5:-5, 5:-5]
        fig = plt.figure()
        contour_local(fig, levels, "\\rho", x, y, rho)
        # plt.savefig("./images/kh_image%s.png" % num, bbox_inches='tight', dpi=300)
        plt.savefig("kh_output.pdf", bbox_inches='tight', dpi=300)
        plt.clf()
        fig.clf()
        f.close()


# fnums = [i*100 for i in range(1,601)]
# fnames = ["./movies/opensbli_%d.h5" % i for i in fnums]
fnums = [0]
# fnames = ["./movies/opensbli_45000.h5"]
fnames = ["./opensbli_output.h5"]
files = zip(fnums, fnames)
plot(files, 20, 0.5, 2.4)
