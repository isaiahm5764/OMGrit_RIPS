array=(64 128 192 256 384 512 768 1024)
for i in "${array[@]}"
	do
	mpirun -np 6 /home/local/IPAMNET/toverman/omgrit-rips/OMGrit_RIPS/examples/advec-diff-imp -mi 10000 -ntime $i -nu 0.1 -mspace 16 -tol 1e-6 -ml 2
	done