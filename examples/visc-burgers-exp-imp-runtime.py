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
time=[200, 256, 512, 1024, 2048, 4096, 5000]
space=[16]
nu=['0.050000']
alpha=['0.500000'] 
ml=[2, 3]
count = 0

exp_ml2=zeros(len(time))
exp_ml3=zeros(len(time))
imp_ml2=zeros(len(time))
imp_ml3=zeros(len(time))
"""
#Open a text file to construct the omgrit matrix in
file_out_conv = open('visc-burgers-newt-conv.txt','w')
file_out_div = open('visc-burgers-newt-div.txt','w')
"""
for t in time:
    for s in space:
        for n in nu:
            for a in alpha:
                for m in ml:
                    with open('out/visc-burgers-newt.time.'+str(t)+"."+str(s)+"."+n+"."+a+"."+str(m)) as f:
                        lines=f.readlines()
                    if(m==2):
                        imp_ml2[count]=lines[0]
                    else:
                        imp_ml3[count]=lines[0]

                    with open('out/visc-burgers-exp-rms.time.'+str(t)+"."+str(s)+"."+n+"."+a+"."+str(m)) as f:
                        lines=f.readlines()
                    if(m==2):
                        exp_ml2[count]=lines[0]
                    else:
                        exp_ml3[count]=lines[0]
    count = count + 1

"""             
file_out_conv.close()
file_out_div.close() 
"""
"""
count = 0
for t in time:
    print(str(imp_ml2[count])+" ")
#print("\n\n")

count = 0
for t in time:
    print(str(imp_ml3[count])+" ")
#print("\n\n")

count = 0
for t in time:
    print(str(exp_ml2[count])+" ")
#print("\n\n")

count = 0
for t in time:
    print(str(exp_ml3[count])+" ")
#print("\n\n")
"""

mpl.figure(1)
mpl.plot(time, exp_ml2, "r")
mpl.plot(time, exp_ml3, "--r")
mpl.plot(time, imp_ml2, "b")
mpl.plot(time, imp_ml3, "--b")
mpl.xlabel('$1/\Delta t$')
mpl.ylabel('Run Time (s)')
mpl.title('Run times for Implicit and Explicit Nonlinear Problem')
mpl.legend(['Exp, 2 level', 'Exp, 3 level',"Imp, 2 level", "Imp, 3 level"])
mpl.show()