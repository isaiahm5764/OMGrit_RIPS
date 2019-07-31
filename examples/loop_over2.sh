array=(1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130)
for i in "${array[@]}"
	do
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-omgrit-full-res -mi $i -ntime 1024 -gamma .5 -tol 0.0
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gs_ex1 -mi $i -ntime 1024 -gamma .5 -tol 0.0
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gj_ex1 -mi $i -ntime 1024 -gamma .5 -tol 0.0
	done