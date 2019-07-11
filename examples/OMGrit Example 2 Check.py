from gurobipy import *
import math

tstart = 0.0
tstop = 1.0
n = 100
deltat = (tstop-tstart)/n

sstart = 0.0
sstop = 1.0
m = 100
deltas = (sstop-sstart)/m

nu = 0.7
alpha = 0.05

A = -(nu*deltat)/(deltas**2) - deltat/(2*deltas)
B = 1 + (2*nu*deltat)/(deltas**2)
C = deltat/(2*deltas) - (nu*deltat)/(deltas**2)

Unaught = [1 for i in range(m/2)] + [0 for i in range(m/2, m+1)]

#KUnaught = [Unaught[1]*B+Unaught[2]*C]
#for i in range(2,m-1):
#	KUnaught.append(Unaught[i-1]*A+Unaught[i]*B+Unaught[i+1]*C)
#KUnaught.append(Unaught[m-2]*A+Unaught[m-1]*B)

try:
	# Initialize model
	M = Model("lp")

	# Import variables
	Uvars = M.addVars(((i,j) for i in range(n+1) for j in range(m+1)), name='U')
	Vvars = M.addVars(((i,j) for i in range(1,n+1) for j in range(1,m)), name='V')

	for j in range(1,m):
		Uvars[0,j].UB = Unaught[j] # initial conditions
		Uvars[0,j].LB = Unaught[j]
	for i in range(n+1):
		Uvars[i,0].UB == 0 # boundary conditions
		Uvars[i,0].LB == 0
		Uvars[i,m].UB == 0
		Uvars[i,m].LB == 0

	# add constraints
	M.addConstr(B*Uvars[1,1] + C*Uvars[1,2] - deltat*Vvars[1,1] == Unaught[1])
	M.addConstrs(A*Uvars[1,j-1] + B*Uvars[1,j] + C*Uvars[1,j+1] - deltat*Vvars[1,j] == Unaught[j] for j in range(2, m-1))
	M.addConstr(A*Uvars[1,m-2] + B*Uvars[1,m-1] - deltat*Vvars[1,m-1] == Unaught[m-1])
	for i in range(2,n+1):
		M.addConstr(-Uvars[i-1,1] + B*Uvars[i,1] + C*Uvars[i,2] - deltat*Vvars[i,1] == 0)
		M.addConstrs(-Uvars[i-1,j] + A*Uvars[i,j-1] + B*Uvars[i,j] + C*Uvars[i,j+1] - deltat*Vvars[i,j] == 0 for j in range(2, m-1))
		M.addConstr(-Uvars[i-1,m-1] + A*Uvars[i,m-2] + B*Uvars[i,m-1] - deltat*Vvars[i,m-1] == 0)

	# Set Objective
	M.setObjective(0.5*deltat*deltas*( sum((Uvars[i,j] - Unaught[j])*(Uvars[i,j] - Unaught[j]) for j in range(m+1) for i in range(n+1)) + alpha*sum(Vvars[i,j]*Vvars[i,j] for j in range(1,m) for i in range(1,n+1))))

	# Prints constraints
	#M.computeIIS()
	M.write("Model.lp")
	
	# Solves model
	#M.Params.Presolve = 0
	M.optimize()
	


	# Writes out solution
	docu = "" # initializes document
	docv = ""
	s = 0 # initialize space
	t = 0 # initialize time
	flag = False
	for v in M.getVars():
		if flag:
			docv += str(v.x) + '\n'
			s += 1
			if s == m:
				s = 1
				t += 1
		else:
			docu += str(v.x) + '\n'
			s += 1
			if s == m+1:
				s = 0
				t += 1
			if t == n+1:
				flag = True
				t = 1
				s = 1

	with open('OMGritEx2SolvedU.txt', 'w') as txtfile:
		txtfile.write(docu)
	with open('OMGritEx2SolvedV.txt', 'w') as txtfile:
		txtfile.write(docv)
except GurobiError:
	print "Gurobi Error"