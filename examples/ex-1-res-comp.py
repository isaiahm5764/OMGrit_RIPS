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

ntime = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125,130];
ex_04_gs = zeros(len(ntime))
ex_04_gj = zeros(len(ntime))
ex_04_omgrit = zeros(len(ntime))

count=0
for n in ntime:
    with open('out/ex-04-omgrit-full-res.res.' + str(n)) as f:
        lines = f.readlines()
    ex_04_omgrit[count] = lines[0]
    with open('out/block_gs_ex1.res.' + str(n)) as f:
        lines = f.readlines()
    ex_04_gs[count] = lines[0]
    with open('out/block_gj_ex1.res.' + str(n)) as f:
        lines = f.readlines()
    ex_04_gj[count] = lines[0]

    count+=1

mpl.figure(1)
mpl.plot(ntime, ex_04_gs, '-k')
mpl.plot(ntime, ex_04_gj, '-g')
mpl.plot(ntime, ex_04_omgrit, '-b')
mpl.xlabel('Number of iterations')
mpl.ylabel('KKT System Residual')
mpl.title('Run times for Linear Quadratic Optimal Control Problem')
mpl.legend(['Block Gauss Seidel', 'Block Gauss Jacobi',"OMGrit"])

mpl.show()