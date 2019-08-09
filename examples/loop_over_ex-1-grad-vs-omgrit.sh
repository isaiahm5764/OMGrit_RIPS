array=(8 16 32 64 128 256 512 1024 2048)
for i in "${array[@]}"
	do
	mpirun -np 8 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-omgrit -mi 10000 -ntime $i -access 1
	#mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04 -mi 10000 -ntime $i
	#mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-serial -maxiter 10000 -ntime $i
	done