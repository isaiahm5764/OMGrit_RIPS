#time=(64 128 256 512 1024 2048 4096)
time=(256)
space=(6 8 10 12 14 16 18 20)
nu=(.2 .4 .6 .8 1.0 1.2 1.4 1.6 1.8 2.0)
alpha=(.005 .5 1.0) 
ml=(2 3 4 5 6)
for t in "${time[@]}"
	do
		for s in "${space[@]}"
		do
			for n in "${nu[@]}"
			do
				for a in "${alpha[@]}"
				do
					for m in "${ml[@]}"
					do
						mpirun -np 1 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/visc-burgers-newt -mi 200 -ntime $t -mspace $s -nu $n -alpha $a -ml %m	
					done
				done
			done
		done
	done