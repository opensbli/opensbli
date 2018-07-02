import numpy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import h5py
import os.path
import matplotlib.cm as cm
import os

plt.style.use('classic')


class plotFunctions(object):
    def __init__(self):
        return

    def contour_local(self, fig, levels0, label, variable):
        ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
        ax1.set_xlabel(r"$x_0$", fontsize=20)
        ax1.set_ylabel(r"$x_1$", fontsize=20)
        CS = ax1.contourf(self.x, self.y, variable, levels=levels0, cmap=cm.jet)
        divider = make_axes_locatable(ax1)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        ticks_at = numpy.linspace(levels0[0], levels0[-1], 10)
        cbar = plt.colorbar(CS, cax=cax1, ticks=ticks_at, format='%.3f')
        cbar.ax.set_ylabel(r"$%s$" % label, fontsize=20)
        return

    def read_file(self, fname):
        f = h5py.File(fname, 'r')
        group = f["opensbliblock00"]
        return f, group

    def read_dataset(self, group, dataset):
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


class Plot(plotFunctions):
    def __init__(self):
        self.Minf = 2.0
        return

    def extract_coordinates(self):
        fname = 'opensbli_output.h5'
        f, group1 = self.read_file(fname)
        x = self.read_dataset(group1, "x0_B0")
        y = self.read_dataset(group1, "x1_B0")
        dx, dy = x[0, 1], y[1, 0]
        print("Grid size (x,y)  is: (%f, %f)" % (x.shape[1], x.shape[0]))
        print("First grid point dx: %f, dy: %f" % (dx, dy))
        dx, dy = x[0, -1]-x[0, -2], y[-1, 0]-y[-2, 0]
        print("Last grid point dx: %f, dy: %f" % (dx, dy))
        return x, y

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhov = self.read_dataset(group, "rhou1_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        u = rhou/rho
        v = rhov/rho
        p = (0.4)*(rhoE - 0.5*(u**2+v**2)*rho)
        return rho, u, v, rhoE, p

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)
        self.x, self.y = self.extract_coordinates()
        rho, u, v, rhoE, p = self.extract_flow_variables(group1)
        variables = [rho, u, v, rhoE, p]
        names = ["rho", "u", "v", "rho E", "P"]

        # Contour plots
        for var, name in zip(variables, names):
            min_val = numpy.min(var)
            max_val = numpy.max(var)
            levels = numpy.linspace(min_val, max_val, n_levels)
            print("%s" % name)
            print(levels)
            fig = plt.figure()
            self.contour_local(fig, levels, "%s" % name, var)
            plt.savefig(directory + "output_%s.pdf" % name, bbox_inches='tight')
            plt.clf()
        f.close()


fname = "opensbli_output.h5"
n_contour_levels = 25
directory = './simulation_plots/'

if not os.path.exists(directory):
    os.makedirs(directory)

KP = Plot()
KP.main_plot(fname, n_contour_levels)
