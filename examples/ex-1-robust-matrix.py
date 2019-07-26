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
gamma=['0.005000', '0.054750', '0.104500', '0.154250', '0.204000', '0.253750', '0.303500', '0.353250', '0.403000', '0.452750', '0.502500', '0.552250', '0.602000', '0.651750', '0.701500', '0.751250', '0.801000', '0.850750', '0.900500', '0.950250', '1.000000']

#Open a text file to construct the omgrit matrix in
file_out_omgrit = open('ex-04-omgrit-robust-matrix.txt','w')
for n in ntime:
    for g in gamma:
        with open('out/ex-04-omgrit-full-res.conv.'+str(n)+"."+g) as f:
            lines=f.readlines()
        file_out_omgrit.write(str(lines[0])+" ")
    file_out_omgrit.write("0;\n")
for g in gamma:
    file_out_omgrit.write("0 ")
file_out_omgrit.close() 

#Open a text file to construct the gs matrix in
file_out_gs = open('ex-04-gs-robust-matrix.txt','w')
for n in ntime:
    for g in gamma:
        with open('out/block_gs_ex1.conv.'+str(n)+"."+g) as f:
            lines=f.readlines()
        file_out_gs.write(str(lines[0])+" ")
    file_out_gs.write("0;\n")
for g in gamma:
    file_out_gs.write("0 ")
file_out_gs.close() 

#Open a text file to construct the gs matrix in
file_out_gj = open('ex-04-gj-robust-matrix.txt','w')
for n in ntime:
    for g in gamma:
        with open('out/block_gj_ex1.conv.'+str(n)+"."+g) as f:
            lines=f.readlines()
        file_out_gj.write(str(lines[0])+" ")
    file_out_gj.write("0;\n")
for g in gamma:
    file_out_gj.write("0 ")
file_out_gj.close() 