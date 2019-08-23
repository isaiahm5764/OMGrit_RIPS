time=(200 256 512 1024 2048 4096 5000)
space=(12)
nu=(.01 .05 .1 .2 .4 .6 .8 1.0 1.2 1.4 1.6 1.8 2.0)
#nu=(.05)
alpha=(.5) 
#ml=(2 3 4 5 6)
ml=(3)
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
						mpirun -np 8 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/visc-burgers-exp-rms -mi 400 -ntime $t -mspace $s -nu $n -alpha $a -ml $m	
						mpirun -np 8 /home/local/IPAMNET/imeyers/OMGrit_RIPS/examples/visc-burgers-newt -mi 400 -ntime $t -mspace $s -nu $n -alpha $a -ml $m
					done
				done
			done
		done
	done