for i in {1..4}
do
	for j in "" "_bc"
	do
		for k in "" "_nointra"
		do
			python loss_novs_v_vs.py rem_reps_vs_stacked${k}${j}_${i}/ rem_reps_novs_stacked${k}${j}_${i}/
			python hist_novs_v_vs.py rem_reps_vs_stacked${k}${j}_${i} rem_reps_novs_stacked${k}${j}_${i}
		done
	done
done
