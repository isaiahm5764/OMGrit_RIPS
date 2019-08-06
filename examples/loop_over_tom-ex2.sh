array=(64 128 256 512 1024)
for i in "${array[@]}"
	do
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-diff-imp -mi 10000 -ntime $i -access 1 -nu 0.1 -mspace 16 -tol 1e-6 -ml 2
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-grad-serial -maxiter 50000 -ntime $i -mspace 16 -stepsize 1.0
	done