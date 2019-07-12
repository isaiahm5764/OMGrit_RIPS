from scipy import *
from matplotlib import pyplot as mpl
from os import sys
#get nsteps for time from the u file

mspace=100
ntime=100
#create the t mesh
tmesh = linspace(0,1.0,ntime+1)
xmesh = linspace(0,1.0,mspace+1)
##

state_vec = empty([ntime+1, mspace+1])
v_vec = empty([ntime+1, mspace+1])

counti=0
countj=0
with open('OMGritEx2SolvedU.txt') as f:
    linesU = f.readlines()
#with open('OMGritEx2SolvedV.txt') as f:
#    linesV = f.readlines()
for i in range(ntime+1):
    for j in range(mspace+1):
        state_vec[i][j] = linesU[mspace*i + j][:]
#for i in range(1,ntime+1):
#    for j in range(1,mspace):
#        v_vec[i][j] = linesV[(mspace-1)*(i-1)+j-1][:]

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data
nx, ny = mspace+1, ntime+1
x = range(nx)
y = range(ny)


hf = plt.figure(1)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, state_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('U value')

hf = plt.figure(2)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, v_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('V value')

plt.show()