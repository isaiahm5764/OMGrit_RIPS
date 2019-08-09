array=(64 128 192 256 384 512 768 1024)
for i in "${array[@]}"
	do
<<<<<<< HEAD
	mpirun -np 1 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-diff-imp -mi 10000 -ntime $i -access 1 -nu 0.1 -mspace 16 -tol 1e-6 -ml 2
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-grad-serial -maxiter 50000 -ntime $i -mspace 16 -stepsize 1.0
=======
	mpirun -np 6 /home/local/IPAMNET/toverman/omgrit-rips/OMGrit_RIPS/examples/advec-diff-imp -mi 10000 -ntime $i -nu 0.1 -mspace 16 -tol 1e-6 -ml 2
>>>>>>> e87f75b6eda89ee3d231f41534aa86cda8308e95
	done