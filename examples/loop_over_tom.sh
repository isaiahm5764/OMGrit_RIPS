array=(8 16 32 64 128 256 512 1024 2048)
for i in "${array[@]}"
	do
	mpirun -np 10 /home/local/IPAMNET/toverman/omgrit-rips/OMGrit_RIPS/examples/ex-04-omgrit -mi 10000 -ntime $i -access 1
	done