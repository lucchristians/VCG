#!/bin/bash
for i in "Qwt" 
do
		file=${i}'_us'
                curr=$(pwd)
                path=${curr}/${file}
                direction=3 #y
		pullgroups=" 6 REST PULLT y RESB PULLB y"
		
		
		#continue
		bash get_distances.sh ${path} ${direction}
		bash setup_us_files.bash ${path}
		for l in {0..29}
		do
			input="${path}/d_${l}/ ${l}${pullgroups}"
			#python setup_us_plummed.py ${input}
		       	
		done
		bash position_r0.bash ${path} 0.0 0.59 0.0 0.0 0.88 0.0
done
bash us_mdp.bash
