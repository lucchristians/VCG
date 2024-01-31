for i in "vs" "novs"
do
	for j in "" "_nointra"
	do
		for k in "" "_bc"
		do
			cp pairhist.txt rem_reps_${i}_stacked${j}${k}/
			cp pairhist_AA.txt rem_reps_${i}_stacked${j}${k}/reference_AA/
			cp const_AA.txt rem_reps_${i}_stacked${j}${k}/reference_AA/
			cp AA_hist.txt rem_reps_${i}_stacked${j}${k}/reference_AA/

			if [ ${i} == "vs" ]
			then
				cp /mnt/data0/REM/setup_vs/vs_r0_stacked/outfile.in${k}${j}.settings rem_reps_${i}_stacked${j}${k}/outfile.in.settings
				cp /mnt/data0/REM/setup_vs/vs_r0_stacked/outfile.data rem_reps_${i}_stacked${j}${k}/
				cp /mnt/data0/REM/setup_vs/vs_r0_stacked/K_const_ind.txt rem_reps_${i}_stacked${j}${k}/
			else
				cp /mnt/data0/REM/setup_vs/novs_stacked/outfile.in${k}${j}.settings rem_reps_${i}_stacked${j}${k}/outfile.in.settings
				cp /mnt/data0/REM/setup_vs/novs_stacked/outfile.data rem_reps_${i}_stacked${j}${k}/
			fi
		done
	done
done
