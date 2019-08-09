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

ntime = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
proc = [1,2,4,6,8,10]
#ntime = [64, 128, 256, 512]
advec_diff_imp = zeros([len(proc),len(ntime)])

pcount=0
for p in proc:
	count=0
	for n in ntime:
	    with open('out/ex-04.time.' + str(n) + '.' + str(p)) as f:
	        lines = f.readlines()
	    advec_diff_imp[pcount,count] = lines[0]
	    count+=1
	pcount+=1


mpl.figure(1)
mpl.plot(ntime, advec_diff_imp[0], '-')
mpl.plot(ntime, advec_diff_imp[1], '--')
mpl.plot(ntime, advec_diff_imp[2], '-.')
mpl.plot(ntime, advec_diff_imp[3], ':')
mpl.plot(ntime, advec_diff_imp[4], '--', lw = 2)
mpl.plot(ntime, advec_diff_imp[5], '-.', lw = 2)
mpl.xlabel('Number of time steps')
mpl.ylabel('Run Time (s)')
mpl.title('Run times for Quadratic Optimal Control Problem ($\\alpha=0.1$)')
mpl.legend(['1 Processor', '2 Processors', '4 Processors', '6 Processors', '8 Processors', '10 Processors'])

mpl.show()