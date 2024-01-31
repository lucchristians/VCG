path=$1
r0x0=$2
r0y0=$3
r0z0=$4
r0x1=$5
r0y1=$6
r0z1=$7
dx=0
dy=0.1
dz=0
for l in {0..29}
do
	cp ${path}/d_${l}/us_test.dat ${path}/d_${l}/us.dat
	for i in "x" "y" "z"
	do
		for j in "0" "1"
		do
			if [ ${i} == "x" ]
			then
				num0=$(bc <<< "scale=3; ${r0x0}+${l}*${dx}")
				num1=$(bc <<< "scale=3; ${r0x1}+${l}*${dx}")
			elif [ ${i} == "y" ]
			then
				num0=$(bc <<< "scale=3; ${r0y0}+${l}*${dy}")
				num1=$(bc <<< "scale=3; ${r0y1}+${l}*${dy}")
			elif [ ${i} == "z" ]
			then
				num0=$(bc <<< "scale=3; ${r0z0}+${l}*${dz}")
				num1=$(bc <<< "scale=3; ${r0z1}+${l}*${dz}")
			fi
			
			if [ ${j} == "0" ]
			then
				sed -i "s/{dist${i}${j}}/${num0}/g" ${path}/d_${l}/us.dat
			elif [ ${j} == "1" ]
			then
				sed -i "s/{dist${i}${j}}/${num1}/g" ${path}/d_${l}/us.dat
			fi
		done
	done
done
