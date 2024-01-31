runs=(rem_reps_*/)
for run in "${runs[@]}"
do
	num=1
	until [ ! -d ${run}rem_iteration${num} ]
	do
		((num++))
	done
	((num--))
	for j in {1..8}
	do
		if [ ! -f ${run}rem_iteration${num}/sub${j}.lammpstrj ]
		then
			echo "${run} ${j} in iter${num} does not exist"
		fi
	done
done
