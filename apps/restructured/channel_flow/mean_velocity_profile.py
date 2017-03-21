from numpy import *
import h5py
import matplotlib.pyplot as plt
import argparse
set_printoptions(threshold='nan')
Nx = 64
Ny = 64
Nz = 64

def plot(path_to_hdf5_file):

    f = h5py.File(path_to_hdf5_file, 'r')
    group = f['opensbliblock00']
    u = group["rhou0_B0"].value/group["rho_B0"].value
    u = u[2:Nx+2, 2:Ny+2, 2:Nz+2] # Remove halo points

    # Compute mean velocity.
    y = linspace(0.0, 2.0, Ny)
    ubar = zeros(Ny)
    for j in range(0, Ny):
        for i in range(0, Nx):
            for k in range(0, Nz):
                ubar[j] += u[k][j][i]
        ubar[j] /= Nx*Ny
        
    plt.plot(y, ubar, "k-", label="OpenSBLI", markerfacecolor='None')

    #plt.axis([0, 2, 0, 30])
    plt.legend(loc="best", numpoints=1)
    plt.xlabel(r"$y$")
    plt.ylabel(r"$\overline{u}^+$")
    plt.savefig("mean_velocity_profile.pdf", bbox_inches='tight')
    
    return
    

if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user.
    parser = argparse.ArgumentParser(prog="mean_velocity_profile")
    parser.add_argument("path", help="The path to the HDF5 file.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
