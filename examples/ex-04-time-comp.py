from scipy import *
from matplotlib import pyplot as mpl
from os import sys

##
# Run like 
#  $ python viz-ex-04.py serial|braid
#
#  Vizualize the files from ex-04.c or ex-04-serial.c
#       - ex-04.out.state      [ output from ex-04-serial ]
#       - ex-04.out.***.###    [ output from ex-04, *** is the step number, and ### is the proc number ]
#       - ex-04.out.design     [ output by both ex-04 and ex-04-serial ]
# 
#  The file ex-04.out.design is a file with Ntime lines, and line k has time
#  step k's design value. 
#
#  The file ex-04.out.state is a file with Ntime lines, and line k has time
#  step k's two solution values.  This is only output from ex-04-serial.
#
#  The files ex-04.out.***.###  This file contains one value, the two solution
#  values at the corresponding step number.  This is only output from ex-04.
#
##

ntime = [8, 16, 32, 64, 128, 256, 512];
ex_04 = zeros(len(ntime))
ex_04_serial = zeros(len(ntime))
ex_04_omgrit = zeros(len(ntime))

count=0
for n in ntime:
    with open('out/ex-04.time.' + str(n)) as f:
        lines = f.readlines()
    ex_04[count] = lines[0]
    with open('out/ex-04-serial.time.' + str(n)) as f:
        lines = f.readlines()
    ex_04_serial[count] = lines[0]
    with open('out/ex-04-omgrit.time.' + str(n)) as f:
        lines = f.readlines()
    ex_04_omgrit[count] = lines[0]

    count+=1
nt=0;
with open('ex-04.out.state') as f:
    lines = f.readlines()
    nt = len(lines)
state_serial = empty([nt, 2])
state_omgrit = empty([nt, 2])
state_par = empty([nt, 2])
with open('ex-04.out.state') as f:
    lines = f.readlines()
count = 0
for line in lines:
    split = line.split("  ")
    state_serial[count][0] = float(split[0].replace(" ", ""))
    state_serial[count][1] = float(split[1].replace(" ", ""))
    count+=1

with open('ex-04.out.u.000') as f:
    lines = f.readlines()
count = 0
for line in lines:
    split = line[6:].split(",")
    state_omgrit[count][0] = float(split[0].replace(" ", ""))
    state_omgrit[count][1] = float(split[1].replace(" ", ""))
    count+=1

tmesh = linspace(0,1.0,nt+1)

for i in range(0,nt):
    si = "%04d" % i
    with open('ex-04.out.'+ si + '.000') as f:
        lines = f.readlines()
    split = lines[0].split(",")
    state_par[i][0] = float(split[0].replace(" ", ""))
    state_par[i][1] = float(split[1].replace(" ", ""))

state_actual = empty([nt, 2])
v_vec = empty(nt)
with open('ex-04.out.v.000') as f:
    lines = f.readlines()
count =0
for line in lines:
    split = line[6:]
    v_vec[count]=split
    count+=1
count=0
dt=1./nt
tmesh2 = linspace(0,1.0,nt)
for t in tmesh2:
    state_actual[count][0] = exp(-t) - 1
    count+=1
count=0
tmesh2 = linspace(0,1.0,nt)
for t in tmesh2:
    state_actual[count][1] = -exp(-t)
    count+=1



mpl.figure(1)
mpl.plot(ntime, ex_04_serial, '-k')
mpl.plot(ntime, ex_04_omgrit, '--')
mpl.xlabel('Number of time steps')
mpl.ylabel('Run Time (s)')
mpl.title('Run times for Linear Quadratic Optimal Control Problem ($\gamma=0.005$)')
mpl.legend(['Gradient descent', 'OMGrit'])
"""
mpl.figure(2)
mpl.plot(tmesh[1:], state_serial[:,0], '-k')
mpl.plot(tmesh[1:], state_omgrit[:,0], '-g')
mpl.plot(tmesh[1:], state_par[:,0], '-r')
mpl.plot(tmesh[1:], state_actual[:,0], '-b')
mpl.xlabel('Time (s)')
mpl.ylabel('State value U_1')
mpl.title('State solution U_1')
mpl.legend(['Serial steepest descent solution', 'OMGrit solution', 'Parallel steepest descent solution'])


mpl.figure(3)
mpl.plot(tmesh[1:], state_serial[:,1], '-k')
mpl.plot(tmesh[1:], state_omgrit[:,1], '-g')
mpl.plot(tmesh[1:], state_par[:,1], '-r')
mpl.plot(tmesh[1:], state_actual[:,1], '-b')
mpl.xlabel('Time (s)')
mpl.ylabel('State value U_2')
mpl.title('State solution U_2')
mpl.legend(['Serial steepest descent solution', 'OMGrit solution', 'Parallel steepest descent solution'])
"""
mpl.show()