job=$1
runs=(*/)
for run in "${runs[@]}"
do
	cd ${run}
	pwd
	python ${job}
	
	cd ..
done
