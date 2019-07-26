#ntime_range=(2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384)
ntime_range=(4096)
#nu_range=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0)
nu_range=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2)
for i in "${ntime_range[@]}"
	do
	for j in "${nu_range[@]}"
		do		
		mpirun /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-diff-imp-full-res -mi 100 -ntime $i -nu $j -mspace 12 -ml 3 -tol 1e-6
		#/home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gs -mi 1000 -ntime $i -mspace 12 -nu $j -tol 1e-6
		#/home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gj -mi 8000 -ntime $i -mspace 12 -nu $j -tol 1e-6		
		done
	done