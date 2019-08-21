#Code below is used to visualize the direct PDE solves for viscous burgers and advection diffusion PDEs
#as solved in advec-diff-serial and visc-burgers-serial

from scipy import *
from matplotlib import pyplot as mpl
from os import sys

###### Helper Functions ######
def retrieve_data(file_stem, var, vec, step):
    """
    Extracts solution data from desired files
    
    Input: var, a string denoting which variable
    data is associated with. vec, a 2D list to
    store data. step, an integer denoting subsection
    of data is being retrieved

    Output: vec, a 2D list containing solution data
    """
    fname = file_stem  + var + ".%03d"%step

    with open(fname) as f:
        lines = f.readlines()
        for line in lines:
            index = int(line[:5]) - 1
            line = line[6:]
            split = line.split(',')
            count = 0
            if var == "u":
                vec[index,0] = left_bound
                vec[index,mspace+1] = right_bound
                count = 1
            for val in split:
                vec[index,count] = float(val)
                count += 1

    return vec


def plot_data(title, fig_num, vec, xmesh, tmesh):
    """
    Visualizes solution data
    Input: title, string containing title of plot.
    fig_num, (fig_num)th figure being generated.
    vec, 2D list containing solution to problem.
    xmesh, spatial mesh to plot solution on.
    tmesh, temporal mesh to plot solution on.
    Output: None
    """
    hf = plt.figure(fig_num)
    ha = hf.add_subplot(111, projection='3d')

    X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
    ha.plot_surface(X, Y, vec)
    ha.set_xlabel('X (position)')
    ha.set_ylabel('t (time)')
    ha.set_zlabel(title)



###### Input ######
# Default input
file_stem = "out/visc-burgers-serial.out."
num_procs = 1
start_t = 0
end_t = 1
start_x = 0
end_x = 1
left_bound = 0
right_bound = 0
nsteps = None # added later if user does not define

# User input
arg_index = 1
while arg_index < len(sys.argv):
    if sys.argv[arg_index] == "-help":
        print " Viz.py plots state, control, and adjoint solutions to time and 1D space dependent optimization problem \n\n"
        print "  -file_stem <file_stem>    : Folder/file location of data (e.g. out/optimize.sol.)\n"
        print "  -nsteps <nsteps>          : Num points in time\n"
        print "  -start_t <start_t>        : Initial time point\n"
        print "  -end_t <end_t>            : Final time point\n"
        print "  -start_x <start_x>        : Initial space point\n"
        print "  -end_x <end_x>            : Final space point\n"
        print "  -l_bound  <l_bound>       : Value of initial space point at all times\n"
        print "  -r_bound <r_bound>        : Value of final space point at all times\n"
        exit()

    elif sys.argv[arg_index] == "-file_stem":
        arg_index += 1
        file_stem = sys.argv[arg_index]
        arg_index += 1

    elif sys.argv[arg_index] == "-nsteps":
        arg_index += 1
        nsteps = int(sys.argv[arg_index])
        arg_index += 1

    elif sys.argv[arg_index] == "-start_t":
        arg_index += 1
        start_t = int(sys.argv[arg_index])
        arg_index += 1

    elif sys.argv[arg_index] == "-end_t":
        arg_index += 1
        end_t = int(sys.argv[arg_index])
        arg_index += 1

    elif sys.argv[arg_index] == "-start_x":
        arg_index += 1
        start_x = int(sys.argv[arg_index])
        arg_index += 1

    elif sys.argv[arg_index] == "-end_x":
        arg_index += 1
        end_x = int(sys.argv[arg_index])
        arg_index += 1

    elif sys.argv[arg_index] == "-l_bound":
        arg_index += 1
        left_bound = int(sys.argv[arg_index])
        arg_index += 1

    elif sys.argv[arg_index] == "-r_bound":
        arg_index += 1
        right_bound = int(sys.argv[arg_index])
        arg_index += 1
    else:
        print "ABORTING: incorrect command line parameter " + str(sys.argv[arg_index])
        exit()


###### Vizualization Code ######
# Get nsteps for time and mspace for space from the u file
with open(file_stem + 'u.000') as f:
        lines = f.readlines()
for line in lines:
    line = line[6:]
    split = line.split(',')
mspace=len(split)
if nsteps == None: # estimate nsteps if user did not define
    nsteps=len(lines)*num_procs


# Create the time and space meshes
tmesh_state = linspace(start_t,end_t,nsteps+1)
xmesh_state = linspace(start_x,end_x,mspace+2)

# Initialize solution vectors
state_vec = zeros([nsteps+1, mspace+2])

# Retrieves solutions of u, v, and w for each processor
state_vec = retrieve_data(file_stem, 'u', state_vec, 0)

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data
nx, ny = mspace, nsteps
x = range(nx)
y = range(ny)

plot_data('U value', 1, state_vec, xmesh_state, tmesh_state)

plt.show()