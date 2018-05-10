from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import math
import h5py

fig = plt.figure()
ax = fig.gca()

#Grid Dimensions
Lx=150
Ly=10

#Stretching factor
by=3.19

#Discretization
dx1=0.15    #Coarser upstream and downstream of the roughness
dx2=0.15    #Finer around the roughness
nx=120/dx1+30/dx2+1     #Number of points in x direction
ny=102                  #Number of points in y direction

#Roughness dimensions
W=6.0       #width (length)
Xr=53.0     #streamwise centre location
k=1         #height
Xini=Xr-W   #Start of the roughness
Xfin=Xr+W   #End of the roughness
i_ini=45/dx1+2/dx2               #First index of the roughness element formula
i_end=45/dx1+(Xfin-45)/dx2    #Last index of the roughness element formula

x_tild=np.arange(-W, W+dx2, dx2)

#Points on x-axis
xcol=np.concatenate((np.arange(0,45,dx1),np.arange(45,75,dx2),np.arange(75,150+dx1,dx1)), axis=1)
x=np.zeros((nx,ny))
for i in range(int(nx)):
      for j in range(ny):   
            x[i][j]=xcol[i]

#First points on y-axis (around the roughness element)
y0=np.zeros(nx)

j_factor=1
c3=0.6565
c4=2.28478
c5=6
c1=c3/2*(1+np.cos(j_factor*math.pi/ny))
c2=c4/2*(1+c5*(j_factor-1)/(ny-1))


y0[i_ini:i_end+1]=-c1*(np.tanh(np.sqrt(x_tild*x_tild)/c2-1)+np.tanh(-np.sqrt(x_tild*x_tild)/c2-1))

#Points on y-axis
eta=np.linspace(0,1,ny) #uniform computational grid

y=np.zeros((nx,ny))
for i in range(int(nx)):
      for j in range(ny):
            y[i][j]=y0[i]+(Ly-y0[i])*np.sinh(by*eta[j])/np.sinh(by)

x=np.transpose(x)
y=np.transpose(y)

#Export mesh to hdf5 file
h5f = h5py.File('smooth_bump.h5', 'w')
h5f.create_dataset('x0', data=x)
h5f.create_dataset('x1', data=y)
h5f.close()

#Plot the grid
plt.plot(x, y, '*')
plt.axis([20, 100, 0, 10])
plt.show()

