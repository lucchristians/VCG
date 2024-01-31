#!/bin/bash

sims=(*_us/)
for sim in "${sims[@]}"
do
	if [ -d ${sim} ]
	then
		cd ${sim}
	else
		echo ${sim} does not exist 
		continue
	fi
	for j in "antiparallel" "parallel" "stacked"
	do
		if [ -d ${j} ]
		then
			cd ${j}
		else
			echo ${sim}${j}/ does not exist
			continue
		fi
		dists=(s_*/)
		for dist in "${dists[@]}"
		do
			if [ -d ${dist} ]
			then 
				cd ${dist}
			else
				echo ${sim}${j}/${dist} does not exist
				continue
			fi
			for l in {0..9}
			do
				if [ ! -f d_${l}/d_${l}.gro ]
				then
					echo "${sim}${j}/${dist}${l} not setup"
				fi
			done
			cd ..
		done
		cd ..
	done
	cd ..
done
