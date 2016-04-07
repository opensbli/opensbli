import argparse
import numpy
from math import pi, exp, cos, sin
import matplotlib.pyplot as plt
import h5py

def plot(path):
    f = open(path, "r")
    
    iteration = []
    ke = []
    enstrophy = []    
    
    for line in f:
        i, j, k, l = line.split(", ")
        iteration.append(int(i))
        ke.append(float(j)/float(l))
        enstrophy.append(float(k)/float(l))
    
    plt.plot(iteration, enstrophy, "-k", label="Enstrophy")
    plt.show()
    return        
        

if(__name__ == "__main__"):
    # Parse the command line arguments provided by the user
    parser = argparse.ArgumentParser(prog="plot")
    parser.add_argument("path", help="The path to the data file.", action="store", type=str)
    args = parser.parse_args()
    
    plot(args.path)
