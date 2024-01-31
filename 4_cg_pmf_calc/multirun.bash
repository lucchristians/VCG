curr=$(pwd)
for i in "novs" "excl" "repl"
do
	for j in "intra" #"inter"
	do
		cd ${curr}/${i}/${j}
		bash ../../run.bash > /dev/null &
	done
done
