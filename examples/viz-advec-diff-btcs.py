# this file is intended to view the evolution of the solution over successive iterations
from scipy import *
from matplotlib import pyplot as mpl
from os import sys
#get nsteps for time from the u file

nsteps=256
mspace=32
#create the t mesh
tmesh = linspace(0,1.0,nsteps)

xmesh = linspace(0,1.0,mspace)

size = len(range(0,26))
print(size)
v_vec = empty([size, nsteps, mspace])
count = 0;
for i in range(0,26):
    for j in range(0,256):
        sj = "%04d" % j
        si = "%04d" % i
        with open('out/advec-diff-btcs.v.out.' + si + "." + sj) as f:
            lines = f.readlines()
        split = lines[0].split(',')
        count2 = 0
        for thing in split:
            v_vec[count,j,count2] = float(split[count2])
            count2+=1
    count+=1


# mpl.figure(1)
# for x in range(0,mspace-1):
#     mpl.plot(tmesh[1:], state_vec[:,x], '-b')
# mpl.ylabel('U')
# mpl.xlabel('time')
# mpl.title('Solution Values')

# mpl.figure(2)
# for x in range(0,mspace-1):
#     mpl.plot(tmesh[1:], w_vec[:,x], '-b')
# mpl.ylabel('W')
# mpl.xlabel('time')
# mpl.title('Solution Values')

# mpl.figure(3)
# for x in range(0,mspace-1):
#     mpl.plot(tmesh[1:], v_vec[:,x], '-b')
# mpl.ylabel('V')
# mpl.xlabel('time')
# mpl.title('Solution Values')
# mpl.show()

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data
nx, ny = mspace, nsteps
x = range(nx)
y = range(ny)

count=0
for i in range(0, 26):
    hf = plt.figure(i+1)
    ha = hf.add_subplot(111, projection='3d')

    X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
    ha.plot_surface(X, Y, v_vec[count])
    ha.set_xlabel('X (position)')
    ha.set_ylabel('t (time)')
    ha.set_zlabel('W value')
    count+=1

plt.show()