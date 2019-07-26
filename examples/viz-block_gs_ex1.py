from scipy import *
from matplotlib import pyplot as mpl
from os import sys

#get nsteps for time from the u file
with open('out/block_gs_ex1.out.u.000') as f:
        lines = f.readlines()
nsteps=len(lines)
#create the t mesh
tmesh = linspace(0,1.0,nsteps+1)
##
current_rank = 0
state_vec = zeros((nsteps,2))
w_vec = zeros((nsteps,2))
v_vec = zeros(nsteps)

with open('out/block_gs_ex1.out.u.000') as f:
    lines = f.readlines()
count = 0
for line in lines:
    line = line[6:]
    split = line.split(',')
    state_vec[count,0] = float(split[0])
    state_vec[count,1] = float(split[1])
    count+=1

with open('out/block_gs_ex1.out.w.000') as f:
    lines = f.readlines()
count = 0
for line in lines:
    line = line[6:]
    split = line.split(',')
    w_vec[count,0] = float(split[0])
    w_vec[count,1] = float(split[1])
    count+=1

with open('out/block_gs_ex1.out.v.000') as f:
    lines = f.readlines()
count = 0
for line in lines:
    line = line[6:]
    v_vec[count] = float(line)
    count+=1


mpl.figure(1)
mpl.plot(tmesh[1:], state_vec[:,0], '-b')
mpl.plot(tmesh[1:], state_vec[:,1], '-k')
mpl.ylabel('u')
mpl.xlabel('time')
mpl.title('Solution Values')
mpl.legend(['Component 1', 'Component 2'])

mpl.figure(2)
mpl.plot(tmesh[1:], w_vec[:,0], '-b')
mpl.plot(tmesh[1:], w_vec[:,1], '-k')
mpl.ylabel('w')
mpl.xlabel('time')
mpl.title('Solution Values')
mpl.legend(['Component 1', 'Component 2'])

mpl.figure(3)
mpl.plot(tmesh[1:], v_vec[:], '-b')
mpl.ylabel('v')
mpl.xlabel('time')
mpl.title('Solution Values')
mpl.show()