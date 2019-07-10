from gurobipy import *
import math

tstart = 0.0
tstop = 1.0
n = 2048
deltat = (tstop-tstart)/n

sstart = 0.0
sstop = 1.0
m = 10
deltas = (sstop-sstart)/m

nu = 2
alpha = 0.005

A = (nu*deltat)/(deltas**2) + deltat/(2*deltas)
B = 1 - (2*nu*deltat)/(deltas**2)
C = (nu*deltat)/(deltas**2) - deltat/(2*deltas)

try:
	# Initialize model
	M = Model("lp")

	# Import variables
	Uvars = M.addVars(((i,j) for i in range(n+1) for j in range(m+1)), name='U')
	Vvars = M.addVars(((i,j) for i in range(n+1) for j in range(m+1)), name='V')

	# add constraints
	M.addConstrs(Uvars[0,j] == 1 for j in range(m+1)) # initial conditions
	M.addConstrs(Uvars[i,0] == 0 for i in range(1,n+1)) # boundary conditions
	M.addConstrs(Uvars[i,m] == 0 for i in range(1,n+1))

	M.addConstr(Uvars[1,1] + Uvars[1,2] - deltat*Vvars[1,1] == B+C)
	M.addConstrs(Uvars[1,j-1] + Uvars[1,j] + Uvars[1,j+1] - deltat*Vvars[1,j] == A+B+C for j in range(2, m))
	M.addConstr(Uvars[1,m-1] + Uvars[1,m] - deltat*Vvars[1,m] == A+B)
	for i in range(2,n+1):
		M.addConstr(-B*Uvars[i-1,1] - C*Uvars[i-1,2] + Uvars[i,1] + Uvars[i,2] - deltat*Vvars[i,1] == 0)
		M.addConstrs(-A*Uvars[i-1,j-1] - B*Uvars[i-1,j] - C*Uvars[i-1,j+1] + Uvars[i,j-1] + Uvars[i,j] + Uvars[i,j+1] - deltat*Vvars[i,j] == 0 for j in range(2, m))
		M.addConstr(-A*Uvars[i-1,m-1] - B*Uvars[i-1,m] + Uvars[i,m-1] + Uvars[i,m] - deltat*Vvars[i,m] == 0)

	# Set Objective
	M.setObjective(0.5*deltat*deltas*sum(Uvars[i,j]*Uvars[i,j] - 2*Uvars[i,j] + 1 + alpha*Vvars[i,j]*Vvars[i,j] for j in range(m+1) for i in range(n+1)))

	# Prints constraints
	M.write("Model.lp")

	# Solves model
	M.optimize()

	# Writes out solution
	docu = "" # initializes document
	docv = ""
	s = 0 # initialize space
	t = 0 # initialize time
	flag = True
	for v in M.getVars():
		if flag:
			docu += str(v.x) + '\n'
		else:
			docv += str(v.x) + '\n'
		s += 1
		if s >= m:
			s = 0
			t += 1
		if t == 2048:
			flag = False
			t = 0

	with open('OMGritEx2SolvedU.txt', 'w') as txtfile:
		txtfile.write(docu)
	with open('OMGritEx2SolvedV.txt', 'w') as txtfile:
		txtfile.write(docv)
except GurobiError:
	print "Gurobi Error"