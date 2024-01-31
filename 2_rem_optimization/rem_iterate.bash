#!/bin/bash

#LMP=/home/lchristians/progs/lammps/src/lmp_serial
#LMP=~/Documents/LFC/lmp_serial
LMP=$1 
#number of iterations
#END=250
#END=1
A_const=$2
B_const=$3
gaussinit=(-3.1038 0.6290 -2.7359 0.4074 -2.5972 0.4599 -2.0324 0.4758 -1.4802 0.5757 ${A_const} ${B_const} ${A_const} ${B_const} ${A_const} ${B_const})
reps=$4
END=$5
# format: A1 B1 A2 B2 A3 B3 ...
#gaussinit=( 9.5 1.0 )
#gaussinit=( 110.0 0.4 )

#gaussinit=( 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4)
#gaussinit=( 1.0500 15.1042 1.0087 0.8500 10.8793 0.2705 2.6791 11.1413 0.0308 2.3064 11.7326 0.7632 2.6351 11.8748 0.7570 0.6500 12.4901 0.9347 0.6500 11.5518 0.7797 2.5534 11.5757 0.7934 2.4577 11.2823 0.7275 0.9500 15.1042 1.0602)

#gaussinit=(0.0140 5.0000  0.4000   0.219  5.0000  0.4000  0.302  5.0000  0.4000  1.140 5.0000  0.4000  1.1560  5.0000  0.4000  0.0285  5.0000  0.4000  0.061 5.0000  0.4000  1.137 5.0000  0.4000  0.3405 5.0000  0.4000  1.231 5.0000  0.4000 )
#gaussinit=(-7.0000  0.6000  -7.0000  0.6000  -7.0000  0.6000 -7.0000  0.6000 -7.0000  0.6000 )
#gaussinit=( 2.0 10.0 0.4 2.0 10.0 0.4 0.391 9.9646 0.2870 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.0 10.0 0.4 2.95 17.49470 1.09110)
append=0
num=1
until [ ! -d rem_iteration${num} ]
do
	((num++))
done
if [ ${num} -gt 1 ] && [ $append -eq 1 ]
then
	#((num++))
	init=${num}
	gaussinit=( $(python rem_iterate_pgd.py) )
else
	init=1
	rm -r rem_iteration*
fi
echo $num
for (( it=${init}; it<=$END; it++ )); do
	mkdir rem_iteration${it}
	cd rem_iteration${it}	
	if [ -f ../input ]
	then
		ln -s ../input .
	fi
	if [ -f ../system.init ]
	then
		ln -s ../system.init . #TODO: commented this out for my system
	fi
	if [ -f ../outfile.in.settings ]
	then
		cp ../outfile.in.settings .
	fi
	if [ -f ../outfile.data ]
	then
		ln -s ../outfile.data .
	fi 

	i=1
	int=1
	printf "iter%d: " ${it}
	for coef in "${gaussinit[@]}" ; do
		if [ $(( int % 2 )) -eq 0 ]; then
			sed -i -z "s/B${i}\n/${coef}\n/g" outfile.in.settings
			printf "%s " "${coef} |"
			((i++))
		#elif [ $(( int % 3 )) -eq 1 ]; then
		#        printf "%s " "${coef}"
		#	sed -i "s/K${i} /${coef} /g" outfile.in.settings
		else
			printf "%s " "${coef}"
			sed -i "s/A${i} /${coef} /g" outfile.in.settings
		fi
		((int++))
	done
	printf "\n"
	#exit
	#TODO: specify path to your version of lammps
	#lmp_serial -i input -l log_it${it}.txt -var SEED ${RANDOM}
	for ((i=1;i<=${reps};i++))
        do
                srun --exclusive -n 1 --mpi=pmi2 $LMP -i input -l log_it${it}_${i}.txt -var SEED ${RANDOM} -var ITER ${i} > /dev/null 2>&1 &
        done
	wait
	cd ..
	gaussinit=( $(python rem_iterate_pgd.py) )
	if [ ${#gaussinit[A]} -eq 0 ]
        then
                echo "Error: failure to obtain constant parameters from rem_iterate.py"
                break
        fi
	
	# If everything works remove the trajectory for non 50 iterations
	if [ ! $(( it % 10 )) -eq 0 ] && [ ! ${it} -eq 1 ]
	then
		rm rem_iteration${it}/sub*.lammpstrj
		#echo ${i}
	fi	
done
