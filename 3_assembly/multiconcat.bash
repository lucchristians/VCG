#!/bin/bash

runs=(*/)
limit=12
num=0
for run in "${runs[@]}"
do
	echo ${run}
	bash concat_files.bash ${run} sizes 
	((num++))
	if [ $num -ge $limit ]
	then
		wait
		num=0
	fi
done
wait

