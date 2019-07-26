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


ntime=[2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
ntimeOmgrit=[16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
nu=['0.000000', '0.100000', '0.200000', '0.300000', '0.400000', '0.500000', '0.600000', '0.700000', '0.800000', '0.900000', '1.000000', '1.100000', '1.200000', '1.300000', '1.400000', '1.500000', '1.600000', '1.700000', '1.800000', '1.900000', '2.000000']

#Open a text file to construct the omgrit matrix in
file_out_omgrit = open('ex-2-omgrit-robust-matrix.txt','w')
for n in ntimeOmgrit:
    for g in nu:
        with open('out/advec-diff-imp-full-res.conv.'+str(n)+"."+g) as f:
            lines=f.readlines()
        file_out_omgrit.write(str(lines[0])+" ")
    file_out_omgrit.write("0;\n")
for g in nu:
    file_out_omgrit.write("0 ")
file_out_omgrit.close() 

"""
#Open a text file to construct the gs matrix in
file_out_gs = open('ex-2-gs-robust-matrix.txt','w')
for n in ntime:
    for g in nu:
        with open('out/block_gs.conv.'+str(n)+"."+g) as f:
            lines=f.readlines()
        file_out_gs.write(str(lines[0])+" ")
    file_out_gs.write("0;\n")
for g in nu:
    file_out_gs.write("0 ")
file_out_gs.close() 


#Open a text file to construct the gs matrix in
file_out_gj = open('ex-2-gj-robust-matrix.txt','w')
for n in ntime:
    for g in nu:
        with open('out/block_gj.conv.'+str(n)+"."+g) as f:
            lines=f.readlines()
        file_out_gj.write(str(lines[0])+" ")
    file_out_gj.write("0;\n")
for g in nu:
    file_out_gj.write("0 ")
file_out_gj.close() 
"""