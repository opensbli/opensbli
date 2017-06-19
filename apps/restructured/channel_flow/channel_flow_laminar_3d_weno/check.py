import numpy
import h5py
f = h5py.File("opensbli.h5", 'r')
group = f['opensbliblock00']
print group["rhou0_B0"].value.shape
print group["rhou0_B0"].value
print numpy.isnan(group["rhou0_B0"].value).any()
