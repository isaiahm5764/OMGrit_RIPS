#Code below is used to visualize the direct PDE solves for viscous burgers and advection diffusion PDEs

from scipy import *
from matplotlib import pyplot as mpl
from os import sys
#get nsteps for time from the u file
with open('visc-berg-seq.out.u.000') as f:
        lines = f.readlines()
for line in lines:
    line = line[6:]
    split = line.split(',')
mspace=len(split)
nsteps=len(lines)
#create the t mesh
tmesh = linspace(0,1.0,nsteps)
xmesh = linspace(0,1.0,mspace+1)
##
current_rank = 0
state_vec = empty([nsteps, mspace+1])

with open('visc-berg-seq.out.u.000') as f:
    lines = f.readlines()
count = 0
for line in lines:
    line = line[6:]
    split = line.split(',')
    count2 = 1
    state_vec[count,0] = 0
    for thing in split:
        state_vec[count,count2] = float(split[count2-1])
        count2+=1
    count+=1


import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data
nx, ny = mspace, nsteps
x = range(nx)
y = range(ny)


hf = plt.figure(1)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, state_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('U value')


plt.show()