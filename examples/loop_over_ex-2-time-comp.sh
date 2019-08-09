array=(1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130)
for i in "${array[@]}"
	do
	#mpirun -np 1 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-diff-imp -mi $i -ntime 8196 -nu 1.5 -mspace 12 -ml 3 -tol 0.0
	mpirun -np 8 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/advec-diff-imp -mi $i -ntime 8196 -nu 1.5 -mspace 12 -ml 3 -tol 0.0 
	#/home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gs -mi $i -ntime 8196 -mspace 12 -nu 1.5 -tol 0.0
	#/home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/block_gj -mi $i -ntime 8196 -mspace 12 -nu 1.5 -tol 0.0
	done