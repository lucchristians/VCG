#!/bin/bash
#ind=$1
runs=(vs_intra_repl_*/ ) #"vs_assemb_pak_lconc/")
#runs=(novs_intra_*/ ) #"vs_assemb_pak_lconc/")
#runs=(vs_intra_repl_*/ novs_intra_*/ ) #"vs_assemb_pak_lconc/")
limit=16
num=0
for ind in "2" #"8" #{8..9}
do
	for run in "vs_intra_repl_2" #"${runs[@]}"
	do
		if [ ! -f ${run}/sub${ind}.lammpstrj ]
		then
			continue	
		fi
		echo ${run} ${ind}
		python rdf_calc.py ${run} sub${ind} ${ind} &
		#python find_cluster_novs.py ${run} sub1 1 > ${run}sizes.txt &
		((num++))
		if [ $num -ge $limit ]
		then
			#wait
			#echo "batch complete"
			num=0
		fi
	done
done
wait

