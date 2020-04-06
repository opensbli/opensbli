import numpy
import matplotlib.pyplot as plt
import h5py
import os.path
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

    def contour_local(self, fig, levels0, label, x, y, variable):
        ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
        ax1.set_xlabel(r"$x$", fontsize=20)
        ax1.set_ylabel(r"$y$", fontsize=20)
        CS = ax1.contour(x, y, variable, levels=levels0, colors='black')
        divider = make_axes_locatable(ax1)
        return


class Plot(plotFunctions):
    def __init__(self):
        self.Lx, self.Ly = 1.0, 0.5
        return

    def line_graphs(self, x, variables, names):
        for i, name in enumerate(names):
            if name is 'Reference':
                plt.plot(x[i], variables[i], label=labels[i], color=colors[i], linestyle='--')
            else:
                plt.plot(x[i], variables[i], label=labels[i], color=colors[i])
        plt.xlabel(r'$x$', fontsize=20)
        plt.ylabel(r'$\%s$' % 'rho', fontsize=20)
        plt.legend(loc='best')
        plt.savefig('viscous_shock_tube_Re%d.pdf' % (Re), bbox_inches='tight')
        plt.clf()
        return

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhov = self.read_dataset(group, "rhou1_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        x = self.read_dataset(group, "x0_B0")
        y = self.read_dataset(group, "x1_B0")
        u = rhou/rho
        v = rhov/rho
        p = (0.4)*(rhoE - 0.5*(u**2 + v**2)*rho)
        return rho, u, v, p, x, y

    def line_compare(self, x, rho):
        # Plot a line along the bottom wall starting from x = 0.3.
        nx = x.shape[1]
        start = 0.3
        x_start = int(start/(self.Lx/float(nx)))
        rho, x = rho[0, x_start:], x[0, x_start:]
        # Data from: "G. Zhou et al. Grid-converged solution and analysis of the unsteady viscous flow in a two-dimensional shock tube. Phys. Fluids 30, 016102 (2018)."
        # (Nx, Ny) = (1500, 750) points used for the reference solution
        data = numpy.loadtxt('reference_Re200.txt')
        x_ref, rho_ref = data[x_start:, 0], data[x_start:, 1]
        x_vars, variables = [x, x_ref], [rho, rho_ref]
        names = ['rho', 'Reference']
        self.line_graphs(x_vars, variables, names)
        return

    def generate_contours(self, x, y, rho, n_levels):
        nx, ny = x.shape[1], y.shape[0]
        start = 0.4
        y_top = 0.25
        x_start, y_start, y_end = int(start/(self.Lx/float(nx))), 0, int(y_top/(self.Ly/float(ny)))
        x, y, rho = x[y_start:y_end, x_start:], y[y_start:y_end, x_start:], rho[y_start:y_end, x_start:]
        fig = plt.figure()
        # contour_levels = numpy.linspace(numpy.min(rho), numpy.max(rho), n_levels)
        contour_levels = numpy.linspace(22.0, 121.0, n_levels)
        self.contour_local(fig, contour_levels, '\rho', x, y, rho)
        plt.savefig('viscous_shock_tube_contours_Re%d.pdf' % Re, bbox_inches='tight')
        return

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)
        rho, u, v, p, x, y = self.extract_flow_variables(group1)
        # Validation plot of the density
        self.line_compare(x, rho)
        # Contour plot of the density
        self.generate_contours(x, y, rho, n_levels)
        return

# default_colourmap = 
labels = ['OpenSBLI', 'Reference']
colors = ['k', 'r']
Re = 200
# directory = "./Re%d/" % Re
fname = "opensbli_output.h5"
PC = Plot()
PC.main_plot(directory + fname, 22)
