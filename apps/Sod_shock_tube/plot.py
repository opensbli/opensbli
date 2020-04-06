import numpy
import matplotlib.pyplot as plt
import h5py
import os.path
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

    def line_graphs(self, x, variables, names):
        for i, name in enumerate(names):
            if name is not 'Exact':
                plt.plot(x[i], variables[i], label=labels[i], color=colors[i], marker='o', linestyle='--', markevery=4)
            else:
                plt.plot(x[i], variables[i], label=labels[i], color=colors[i])
        plt.xlabel(r'$x$', fontsize=20)
        plt.ylabel(r'$\%s$' % 'rho', fontsize=20)
        plt.legend(loc='best')
        plt.savefig('Sod_shock_tube.pdf', bbox_inches='tight')
        plt.clf()
        return

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        x = self.read_dataset(group, "x0_B0")
        u = rhou/rho
        p = (0.4)*(rhoE - 0.5*(u**2)*rho)
        return rho, u, p, x

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)
        rho, u, p, x = self.extract_flow_variables(group1)
        data = numpy.loadtxt('reference.txt')
        x_ref, rho_ref = data[:, 0], data[:, 1]
        x_vars, variables = [x, x_ref], [rho, rho_ref]
        names = ['rho', 'Exact']
        self.line_graphs(x_vars, variables, names)
        return

labels = ['OpenSBLI', 'Exact']
colors = ['r', 'k']

directory = './'
fname = "opensbli_output.h5"
PC = Plot()
PC.main_plot(directory + fname, 25)
