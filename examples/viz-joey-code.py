from scipy import *
from matplotlib import pyplot as mpl
from os import sys
#get nsteps for time from the u file

mspace=10
ntime=2048
#create the t mesh
tmesh = linspace(0,1.0,ntime)
xmesh = linspace(0,1.0,mspace)
##

state_vec = empty([ntime, mspace])
v_vec = empty([ntime, mspace])

counti=0
countj=0
subnum=9
with open('OMGritEx2Solved.txt') as f:
    lines = f.readlines()
for i in range(ntime):
    for j in range(mspace):
        state_vec[i][j] = lines[(i+1) * (j+1)][subnum:]
        v_vec[i][j] = lines[2048*10 + ((i+1) * (j+1))][subnum:]
    if i==10:
        subnum +=1
    if i==100:
        subnum +=1
    if i==1000:
        subnum+=1


import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data
nx, ny = mspace, ntime
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
ha.set_zlabel('U value')

plt.show()