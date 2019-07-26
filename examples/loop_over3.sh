ntime_range=(2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384)
gamma_range=(.005 0.05475 0.1045 0.15425 0.204 0.25375 0.3035 0.35325 0.403 0.45275 0.5025 0.55225 0.602 0.65175 0.7015 0.75125 0.801 0.85075 0.9005 0.95025 1.0)
for i in "${ntime_range[@]}"
	do
	for j in "${gamma_range[@]}"
		do	
		mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/ex-04-omgrit-full-res -mi 10000 -ntime $i -gamma $j -tol 1e-6
		mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gs_ex1 -mi 10000 -ntime $i -gamma $j -tol 1e-6
		mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gj_ex1 -mi 10000 -ntime $i -gamma $j -tol 1e-6
		done
	done