from scipy import *
from matplotlib import pyplot as mpl
from os import sys


ntime = [128, 128, 256, 512, 1024]
#ntime = [64, 128, 256, 512]
advec_diff_imp_1p = zeros(len(ntime))
advec_diff_imp_8p = zeros(len(ntime))
advec_grad_serial = zeros(len(ntime))

count=0
for n in ntime:
    with open('out/advec-diff-imp.time.' + str(n)+".1"+".10") as f:
        lines = f.readlines()
    advec_diff_imp_1p[count] = lines[0]
    with open('out/advec-diff-imp.time.' + str(n)+".8"+".10") as f:
        lines = f.readlines()
    advec_diff_imp_8p[count] = lines[0]    
    with open('out/advec-grad-serial.time.' + str(n)+".10") as f:
        lines = f.readlines()
    advec_grad_serial[count] = lines[0]

    count+=1


mpl.figure(1)
mpl.plot(ntime, advec_diff_imp_1p, '--')
mpl.plot(ntime, advec_diff_imp_8p, '-k')
mpl.plot(ntime, advec_grad_serial, '-.')
mpl.xlabel('Number of time steps')
mpl.ylabel('Run Time (s)')
mpl.title('Run times for Advection-Diffusion Model Problem $\\nu=0.05$, 10 space points')
mpl.legend(['OMGrit 1 proc', 'OMGrit 8 proc', 'Gradient Descent'])

mpl.show()