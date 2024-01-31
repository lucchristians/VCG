for i in "novs" "vs" 
do
	for j in "inter" "intra"
	do
		python cg_us_pmf_plots.py ${i}/${j}/ log_${i} ${i}_${j}_potential.png 0 30
	done
done
