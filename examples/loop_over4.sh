array=(2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384)
for i in "${array[@]}"
	do
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-omgrit-full-res -mi 10000 -ntime $i -gamma .5 -tol 1e-6
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gs_ex1 -mi 10000 -ntime $i -gamma .5 -tol 1e-6
	mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gj_ex1 -mi 10000 -ntime $i -gamma .5 -tol 1e-6
	done