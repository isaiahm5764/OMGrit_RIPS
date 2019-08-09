array=(8 16 32 64 128 256 512 1024 2048)
for i in "${array[@]}"
	do
<<<<<<< HEAD
	mpirun -np 1 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-omgrit -mi 10000 -ntime $i -access 1
=======
	mpirun -np 10 /home/local/IPAMNET/toverman/omgrit-rips/OMGrit_RIPS/examples/ex-04-omgrit -mi 10000 -ntime $i -access 1
>>>>>>> e87f75b6eda89ee3d231f41534aa86cda8308e95
	#mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04 -mi 10000 -ntime $i
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-serial -maxiter 10000 -ntime $i
	done