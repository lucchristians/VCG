runs=(*/)
for run in "${runs[@]}"
do
	#cd ${run}
	#pwd
	python const_plot_vs.py ${run}
	
	#cd ..
done
