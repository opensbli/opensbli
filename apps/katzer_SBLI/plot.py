import numpy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import h5py
import os.path
import matplotlib.cm as cm


plt.style.use('classic')

# Matplotlib settings for publication-ready figures
try:
    f = open(os.path.expanduser('./rcparams.py'), 'r')
    exec(f.read())
except:
    pass


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


class KatzerPlot(plotFunctions):
    def __init__(self):
        self.Minf = 2.0
        self.Re = 950
        self.RefT = 288.0
        self.SuthT = 110.4
        self.Ly = 115.0
        self.Lx = 400.0
        self.scale = 2.31669259
        self.D11 = self.extract_metrics()
        return

    def load_reference_data(self):
        reference = numpy.loadtxt('./reference_data.txt')
        x, cf, p = reference[:, 0], reference[:, 1], reference[:, 2]
        return x, cf, p

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

    def extract_metrics(self):
        fname = 'opensbli_output.h5'
        f, group1 = self.read_file(fname)
        D11 = self.read_dataset(group1, "D11_B0")
        return D11

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhov = self.read_dataset(group, "rhou1_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        u = rhou/rho
        v = rhov/rho
        p = (0.4)*(rhoE - 0.5*(u**2+v**2)*rho)
        M = numpy.sqrt(u**2 + v**2)
        T = 1.4*(self.Minf**2)*p/rho
        mu = self.compute_viscosity(T)
        return rho, u, v, rhoE, p, T, M, mu

    def compute_wall_derivative(self, variable):
        ny = numpy.size(self.y[:, 0])
        Ly = self.Ly
        delta = Ly/(ny-1.0)
        D11 = self.D11[0:6, :]
        var = variable[0:6, :]
        coeffs = numpy.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15])
        coeffs = coeffs.reshape([6, 1])
        dudy = sum(D11*var*coeffs)/delta
        return dudy

    def compute_viscosity(self, T):
        mu = (T**(1.5)*(1.0+self.SuthT/self.RefT)/(T+self.SuthT/self.RefT))
        return mu

    def compute_skin_friction(self, u, mu):
        # Wall viscosity all x points
        mu_wall = mu[0, :]
        dudy = self.compute_wall_derivative(u)
        tau_wall = dudy*mu_wall
        Cf = tau_wall/(0.5*self.Re)
        return Cf

    def SBLI_comparison(self, Cf, P):
        # Skin friction plot
        x, ref_cf, ref_p = self.load_reference_data()

        Rex0 = 0.5*(self.Re/self.scale)**2
        x0 = 0.5*self.Re/self.scale**2
        delta = 0.1
        # Calculate local Reynolds number for all x
        Rex = Rex0 + self.Re*self.x[1, :]

        # plt.plot(markus_x, exact_cf, label='Flat plate exact')
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, ref_cf, color='r', linestyle='--', marker='o', markevery=15, markersize=5, label='Reference')

        ax.plot(self.x[1, :], Cf, color='k', label='OpenSBLI')
        ax.axhline(y=0.0, linestyle='--', color='k')

        ax.set_xlabel(r'$x_0$', fontsize=20)
        ax.set_ylabel(r'$C_f$', fontsize=20)
        # ax.set_title('Skin friction')
        plt.legend(loc="best")
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        fig.savefig('skin_friction.pdf', bbox_inches='tight')
        fig.clf()

        plt.plot(x, ref_p, color='r', linestyle='--', marker='o', markevery=15, markersize=5, label='Reference')
        plt.plot(self.x[1, :], P[0, :]/P[0, 0], color='k', label="OpenSBLI")
        # linestyle='', marker='o',markevery=10)
        plt.xlabel(r'$x_0$', fontsize=20)
        plt.ylabel(r'$\frac{P_w}{P_1}$', fontsize=22)
        plt.title('Normalized wall pressure')
        plt.legend(loc="best")
        plt.savefig("wall_pressure.pdf", bbox_inches='tight')
        plt.clf()
        return

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)

        self.x, self.y = self.extract_coordinates()
        rho, u, v, rhoE, p, T, M, mu = self.extract_flow_variables(group1)
        variables = [rho, u, v, rhoE, p, M, T, mu]
        names = ["rho", "u", "v", "rho E", "P", "M", "T", "mu"]

        # Contour plots
        for var, name in zip(variables, names):
            min_val = numpy.min(var)
            max_val = numpy.max(var)
            levels = numpy.linspace(min_val, max_val, n_levels)
            print("%s" % name)
            print(levels)
            fig = plt.figure()
            self.contour_local(fig, levels, "%s" % name, var)
            plt.savefig("katzer_%s.pdf" % name, bbox_inches='tight')
            plt.clf()
        # Compare to SBLI
        Cf = self.compute_skin_friction(u, mu)
        self.SBLI_comparison(Cf, p)
        f.close()


fname = "opensbli_output.h5"
n_contour_levels = 25


KP = KatzerPlot()
KP.main_plot(fname, n_contour_levels)
