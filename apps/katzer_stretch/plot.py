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
    # p = group["p_B0"].value

    rho = rho[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    grid_points = [np[0] - 2*k, np[1] - 2*k]
    rhou = rhou[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    rhov = rhov[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    x = x[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    y = y[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    rhoE = rhoE[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    # p = p[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]

    u = rhou/rho
    v = rhov/rho
    p = (0.4)*(rhoE - 0.5*(u**2+v**2)*rho)
    a = numpy.sqrt(1.4*p/rho)
    M = u/a
    T = 1.4*4*p/rho
    # Cf = group["Cf_B0"].value[0]
    return x, y, rho, u, v, rhoE, p, M, T


def load_SBLI_data():
    markus_cf = numpy.loadtxt('./SBLI_DATA/markus_cf.txt')
    exact_cf = numpy.loadtxt('./SBLI_DATA/markus_cf_exact.txt')
    markus_x = numpy.loadtxt('./SBLI_DATA/markus_x.txt')
    markus_pwall = numpy.loadtxt('./SBLI_DATA/markus_pwall.txt')
    return markus_cf, exact_cf, markus_x, markus_pwall


def initial_profile(u, x, y, lhalo, rhalo):
    f, group1 = read_file('initial.h5')
    rho_in = group1["rho_B0"].value
    rhou_in = group1["rhou0_B0"].value
    np = group1["rho_B0"].shape
    rho_in = rho_in[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    rhou_in = rhou_in[lhalo:np[0]-rhalo, lhalo:np[1]-rhalo]
    u_in = rhou_in/rho_in

    fact = 50
    bl_end = 100
    for i in range(0, 10):

        plt.plot(u[0:bl_end, i*fact], y[0:bl_end, i*fact], '-k', label="Solution")
        plt.plot(u_in[0:bl_end, i*fact], y[0:bl_end, i*fact], '--r', label="Initial condition")

        plt.title(r'Boundary layer profile at $x_0$ = %d' % x[0, i*fact], fontsize=20)
        plt.xlim([0, 1.2])
        plt.legend(loc="best")
        plt.savefig("bl_velocity_%d.pdf" % i, bbox_inches='tight')
        plt.clf()
    plt.ioff()
    fig, axs = plt.subplots(nrows=3, ncols=2, sharey=True, figsize=(12, 16))
    # axs[0].set_ylabel(r'$x_1$', fontsize=20)
    print u[1, :]
    axs = numpy.ndarray.flatten(axs)

    for i in range(0, 6):
        xval = x[0, (i+1)*fact]
        print xval
        print u[1, (i+1)*fact]
        axs[i].plot(u[0:bl_end, (i+1)*fact], y[0:bl_end, (i+1)*fact], '-k', label="Solution")
        axs[i].plot(u_in[0:bl_end, (i+1)*fact], y[0:bl_end, (i+1)*fact], '--r', label="Initial condition")
        axs[i].xaxis.set_ticks(numpy.linspace(0, 1.2, 4))
        axs[i].set_xlabel(r'$u_0$', fontsize=24)
        axs[i].set_ylabel(r'$x_1$', fontsize=24)
        axs[i].set_title(r'$x_0 = %d$' % xval, fontsize=22)
        for tick in axs[i].xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in axs[i].yaxis.get_major_ticks():
            tick.label.set_fontsize(14)

        # axs[i].set(aspect='equal')
    plt.subplots_adjust(wspace=0.2, hspace=0.4)

    # plt.xlabel(r'$u_0$ velocity', fontsize=20)
    # plt.ylabel(r'$x_1$', fontsize=20)
    fig.savefig("profile_comparison.pdf", bbox_inches='tight')
    f.close()
    return


def line_graphs(x, variable, name):
    if name == "u":
        plt.axhline(y=0.0, linestyle='--', color='k')

    plt.plot(x[1, :], variable)
    plt.xlabel(r'$x_0$', fontsize=20)
    plt.ylabel(r'$%s$ at wall' % name, fontsize=20)
    plt.savefig("wall_%s.pdf" % name, bbox_inches='tight')
    plt.clf()
    return


def SBLI_comparison(x, Cf, P):
    # Skin friction plot
    markus_cf, exact_cf, markus_x, markus_pwall = load_SBLI_data()

    scale = 2.31669259
    Re = 950
    Rex0 = 0.5*(Re/scale)**2
    x0 = 0.5*Re/scale**2
    delta = 0.1
    # Calculate local Reynolds number for all x
    Rex = Rex0 + Re*x[1, :]

    # plt.plot(markus_x, exact_cf, label='Flat plate exact')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(markus_x, markus_cf, color='r', linestyle='--', marker='o', markevery=15, markersize=5,
            label='Reference')

    ax.plot(x[1, :], Cf, color='k', label='OpenSBLI')
    ax.axhline(y=0.0, linestyle='--', color='k')

    ax.set_xlabel(r'$x_0$', fontsize=20)
    ax.set_ylabel(r'$C_f$', fontsize=20)
    # ax.set_title('Skin friction')
    plt.legend(loc="best")
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    fig.savefig('skin_friction.pdf', bbox_inches='tight')
    fig.clf()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(markus_x, markus_cf, label='Reference')
    dx = -0.040
    dy = 0
    offset = transforms.ScaledTranslation(dx, dy,
                                          fig.dpi_scale_trans)
    shadow_transform = ax.transData + offset
    ax.plot(x[1, :], Cf, transform=shadow_transform, lw=1, label='Shifted')
    ax.axhline(y=0.0, linestyle='--', color='k')

    ax.set_xlabel(r'$x_0$', fontsize=20)
    ax.set_ylabel(r'$C_f$', fontsize=20)
    ax.set_title('Shifted skin friction')
    plt.legend(loc="best")
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    fig.savefig('shifted_skin_friction.pdf', bbox_inches='tight')
    fig.clf()

    plt.plot(markus_x, markus_pwall, color='r', linestyle='--', marker='o', markevery=15, markersize=5,
             label='Reference')
    plt.plot(x[1, :], P[0, :]/P[0, 0], color='k', label="OpenSBLI")
    # linestyle='', marker='o',markevery=10)
    plt.xlabel(r'$x_0$', fontsize=20)
    plt.ylabel(r'$\frac{P_w}{P_1}$', fontsize=22)
    plt.title('Normalized wall pressure')
    plt.legend(loc="best")
    plt.savefig("wall_pressure.pdf", bbox_inches='tight')
    plt.clf()

    # Write out Cf and P to text:
    # numpy.savetxt('cf.txt', Cf)
    # numpy.savetxt('x.txt', x[1,:])
    # numpy.savetxt('pwall.txt', P[0,:]/P[0,0])
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
        plt.savefig("katzer_%s.pdf" % name, bbox_inches='tight')
        plt.clf()

    # Compare initial profile
    initial_profile(u, x, y, 5, 5)
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

    # Compare to SBLI
    # SBLI_comparison(x, Cf, P)
    f.close()


fname = "opensbli_output.h5"
plot(fname, 25)

# if(__name__ == "__main__"):
#     plot("./")
