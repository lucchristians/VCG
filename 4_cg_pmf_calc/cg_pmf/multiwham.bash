curr_dir=$(pwd)
for i in "novs" "repl" "excl"
do
	cd ${curr_dir}/${i}
	bash ../prep_wham.bash y
	bash ../run_wham.bash
done
