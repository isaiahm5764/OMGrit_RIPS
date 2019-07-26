# this file is intended to view the evolution of the solution over successive iterations
from scipy import *
from matplotlib import pyplot as mpl
from os import sys
#get nsteps for time from the u file

nsteps=256
mspace=2
#create the t mesh
tmesh = linspace(0,1.0,nsteps+1)

xmesh = linspace(0,1.0,mspace)

size = len(range(0,15))
print(size)
v_vec = empty([size, nsteps, mspace])
count = 0;
for i in range(0,15):
    for j in range(0,256):
        sj = "%04d" % j
        si = "%02d" % i
        with open('out/ex-04-omgrit-iter.out.' + si + "." + sj+".000") as f:
            lines = f.readlines()
        split = lines[0].split(',')
        count2 = 0
        for thing in split:
            v_vec[count,j,count2] = float(split[count2])
            count2+=1
    count+=1


import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data

count=0
for i in range(0, 15):
    mpl.figure(i)
    mpl.plot(tmesh[1:], v_vec[i,:,0], '-b')
    mpl.plot(tmesh[1:], v_vec[i,:,1], '-k')
    mpl.ylabel('W')
    mpl.xlabel('time')
    mpl.title('Adjoint Solution Values')
    mpl.legend(['$w_{1}$', '$w_{2}$'])
    mpl.tick_params(left= False, labelleft = False)
    count+=1

plt.show()