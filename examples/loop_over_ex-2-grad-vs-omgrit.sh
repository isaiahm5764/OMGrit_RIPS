array=(128 256 512 1024)
#array=(1024)
for i in "${array[@]}"
	do
	mpirun -np 1 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-diff-imp -mi 10000 -ntime $i -access 1 -nu 0.05 -mspace 20 -tol 1e-6 -ml 3
	#mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-grad-serial -maxiter 45000 -ntime $i -mspace 20 -stepsize 1.0 -nu .05
	done