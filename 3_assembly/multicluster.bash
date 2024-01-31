#!/bin/bash
#ind=$1
runs=(vs_intra_repl_*/ ) #"vs_assemb_pak_lconc/")
limit=16
num=0
for ind in {1..10}
do
	for run in "${runs[@]}"
	do
		if [ ! -f ${run}/sub${ind}.lammpstrj ]
		then
			continue	
		fi
		echo ${run} ${ind}
		python find_cluster.py ${run} sub${ind} ${ind} > ${run}sizes${ind}.txt &
		#python find_cluster_novs.py ${run} sub1 1 > ${run}sizes.txt &
		((num++))
		if [ $num -ge $limit ]
		then
			wait
			num=0
		fi
	done
done
wait

