curr_dir=$(pwd)
for i in "novs" "repl" "excl"
do
	cd ${curr_dir}/${i}
	bash ../prep_windows.bash 
	bash ../run_windows.bash
done
