from scipy import *
from matplotlib import pyplot as mpl
from os import sys

"""
time=[64, 128, 256, 512, 1024, 2048, 4096]
space=[12]
nu=['0.010000', '0.050000', '0.100000', '0.200000', '0.400000', '0.600000', '0.800000', '1.000000', '1.200000', '1.400000', '1.600000', '1.800000', '2.000000']
alpha=['0.005000', '0.500000', '1.000000'] 
ml=[2, 3, 4, 5, 6]
"""
time=[64, 128, 256, 512, 1024, 2048, 4096]
space=[12, 14, 16, 18, 20, 22, 26, 30]
nu=['0.050000']
alpha=['0.500000'] 
ml=[2, 3]

#Open a text file to construct the omgrit matrix in
file_out_conv = open('visc-burgers-exp-rms-conv.txt','w')
file_out_div = open('visc-burgers-exp-rms-div.txt','w')
for t in time:
    for s in space:
        for n in nu:
            for a in alpha:
                for m in ml:
                    with open('out/visc-burgers-exp-rms.conv.'+str(t)+"."+str(s)+"."+n+"."+a+"."+str(m)) as f:
                        lines=f.readlines()
                    if(lines[0]=="1.000000"):
                        file_out_conv.write("ntime, mspace, nu, alpha, max levels, conv/div:"+str(t)+", "+str(s)+", "+n+", "+a+", "+str(m)+", ")
                        file_out_conv.write(str(lines[0])+"\n\n")
                    else:
                        file_out_div.write("ntime, mspace, nu, alpha, max levels, conv/div:"+str(t)+", "+str(s)+", "+n+", "+a+", "+str(m)+", ")
                        file_out_div.write(str(lines[0])+"\n\n")                       
file_out_conv.close()
file_out_div.close() 