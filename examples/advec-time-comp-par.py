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

proc = [1,2,4,6,8,10]
#ntime = [64, 128, 256, 512]
advec_diff_imp = zeros(len(proc))

pcount=0
for p in proc:
	count=0
	with open('out/advec-diff-imp.time.' + str(2048) + '.' + str(p)) as f:
		lines = f.readlines()
	advec_diff_imp[pcount] = lines[0]
	pcount+=1


mpl.figure(1)
mpl.plot(proc, advec_diff_imp, '-')

mpl.xlabel('Number of processors')
mpl.ylabel('Run Time (s)')
mpl.title('Run times for Advection-Diffusion Model Problem ($\\nu=0.1$)\n2048 time points, 40 space points')

mpl.show()