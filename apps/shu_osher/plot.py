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

    def read_file(self, fname):
        f = h5py.File(fname, 'r')
        group = f["opensbliblock00"]
        return f, group

    def read_dataset(self, group, dataset):
        d_m = group["%s" % (dataset)].attrs['d_m']
        size = group["%s" % (dataset)].shape
        read_start = [abs(d) for d in d_m]
        read_end = [s-abs(d) for d, s in zip(d_m, size)]
        if len(read_end) == 1:
            read_data = group["%s" % (dataset)].value[read_start[0]:read_end[0]]
        elif len(read_end) == 2:
            read_data = group["%s" % (dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1]]
        elif len(read_end) == 3:
            read_data = group["%s" % (dataset)].value[read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
        else:
            raise NotImplementedError("")
        return read_data


class Plot(plotFunctions):
    def __init__(self):
        return

    def line_graphs(self, x, variable, name):
        plt.plot(x, variable)
        plt.xlabel(r'$x_0$', fontsize=20)
        plt.ylabel(r'$%s$' % name, fontsize=20)
        plt.savefig(directory + "output_%s.pdf" % name, bbox_inches='tight')
        plt.clf()
        return

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        u = rhou/rho
        p = (0.4)*(rhoE - 0.5*(u**2)*rho)
        return rho, u, rhoE, p

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)
        rho, u, rhoE, p = self.extract_flow_variables(group1)
        variables = [rho, u, p]
        names = ["rho", "u", "P"]
        x = numpy.linspace(0, 10, rho.size)

        for var, name in zip(variables, names):
            self.line_graphs(x, var, name)
            f.close()


fname = "opensbli_output.h5"
n_contour_levels = 25
directory = './simulation_plots/'

if not os.path.exists(directory):
    os.makedirs(directory)

KP = Plot()
KP.main_plot(fname, n_contour_levels)
