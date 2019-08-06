from scipy import *
from matplotlib import pyplot as mpl
from os import sys
#get nsteps for time from the u file
with open('out/advec-diff-imp.out.u.000') as f:
        lines = f.readlines()
for line in lines:
    line = line[6:]
    split = line.split(',')
mspace=len(split)
nsteps=len(lines)*8+1 +8
#create the t mesh
tmesh = linspace(0,1.0,nsteps)
tmesh_state = linspace(0,1.0,nsteps+1)
xmesh = linspace(0,1.0,mspace)
xmesh_state = linspace(0,1.0,mspace+2)
##
current_rank = 0

file_stem = "out/advec-diff-imp.out."
state_vec = empty([nsteps+1, mspace+2])
w_vec = empty([nsteps, mspace])
v_vec = empty([nsteps, mspace])
# We don't know the MPI ranks ahead of time that generated the files, so we guess :-)
with open('out/advec-diff-imp.out.u0.000') as f:
    lines = f.readlines()
split = lines[0].split(',')
count2=1
state_vec[0,0] = 0
state_vec[0,mspace+1] = 0
for thing in split:
    state_vec[0,count2]=float(split[count2-1])
    count2+=1
for step in range(0,8):
    fname = file_stem  + "u." + "%03d"%step
    fname1 = file_stem + "w."  + "%03d"%step
    fname2 = file_stem + "v."  + "%03d"%step
    with open(fname) as f:
        lines = f.readlines()
    for line in lines:
        index = int(line[:5])
        line = line[6:]
        split = line.split(',')
        state_vec[index,0] = 0
        state_vec[index,mspace+1] = 0
        count2=1
        for thing in split:
            state_vec[index,count2] = float(split[count2-1])
            count2 +=1
    with open(fname1) as f:
        lines = f.readlines()
    for line in lines:
        index = int(line[:5])
        line = line[6:]
        split = line.split(',')
        count2=0
        for thing in split:
            w_vec[index,count2] = float(split[count2])
            count2 +=1
    with open(fname2) as f:
        lines = f.readlines()
    for line in lines:
        index = int(line[:5])
        line = line[6:]
        split = line.split(',')
        count2=0
        for thing in split:
            v_vec[index,count2] = float(split[count2])
            count2 +=1

# with open('out/advec-diff-imp.out.u0.000') as f:
#     lines = f.readlines()
# split = lines[0].split(',')
# count2=1
# state_vec[0,0] = 0
# state_vec[0,mspace+1] = 0
# for thing in split:
#     state_vec[0,count2]=float(split[count2-1])
#     count2+=1

# with open('out/advec-diff-imp.out.u.000') as f:
#     lines = f.readlines()
# count = 1
# for line in lines:
#     line = line[6:]
#     split = line.split(',')
#     state_vec[count,0] = 0
#     state_vec[count,mspace+1] = 0
#     count2 = 1
#     for thing in split:
#         state_vec[count,count2] = float(split[count2-1])
#         count2+=1
#     count+=1

# print(state_vec)
# with open('out/advec-diff-imp.out.w.000') as f:
#     lines = f.readlines()
# count = 0
# for line in lines:
#     line = line[6:]
#     split = line.split(',')
#     count2 = 0
#     for thing in split:
#         w_vec[count,count2] = float(split[count2])
#         count2+=1
#     count+=1

# with open('out/advec-diff-imp.out.v.000') as f:
#     lines = f.readlines()
# count = 0
# for line in lines:
#     line = line[6:]
#     split = line.split(',')
#     count2 = 0
#     for thing in split:
#         v_vec[count,count2] = float(split[count2])
#         count2+=1
#     count+=1

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


hf = plt.figure(1)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh_state, tmesh_state)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, state_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('U value')

hf = plt.figure(2)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, w_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('W value')

hf = plt.figure(3)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, v_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('V value')

plt.show()