from scipy import *
from matplotlib import pyplot as mpl
from os import sys


ntime = [8, 16, 32, 64, 128, 256, 512, 1024, 2048]
#ntime = [64, 128, 256, 512]
omgrit_1p = zeros(len(ntime))
omgrit_8p = zeros(len(ntime))
ex4_grad_serial = zeros(len(ntime))

count=0
for n in ntime:
    with open('out/ex-04-omgrit-time.' + str(n)+".1") as f:
        lines = f.readlines()
    omgrit_1p[count] = lines[0]
    with open('out/ex-04-omgrit-time.' + str(n)+".8") as f:
        lines = f.readlines()
    omgrit_8p[count] = lines[0]    
    with open('out/ex-04-serial.time.' + str(n)) as f:
        lines = f.readlines()
    ex4_grad_serial[count] = lines[0]

    count+=1


mpl.figure(1)
mpl.plot(ntime, omgrit_1p, '--')
mpl.plot(ntime, omgrit_8p, '-k')
mpl.plot(ntime, ex4_grad_serial, '-.')
mpl.xlabel('Number of time steps')
mpl.ylabel('Run Time (s)')
mpl.title('Run times for Linear Quadratic Optimal Control Problem $\\gamma=0.005$')
mpl.legend(['OMGrit 1 proc', 'OMGrit 8 proc', 'Gradient Descent'])

mpl.show()