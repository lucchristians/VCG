#!/bin/bash
#Variables

for i in "Qwt" "Q5"
do
        for j in  "s104" "p38" "p28" "a18" "a48"
        do
                file=${i}'_us'
                s=${j:0:1}
                dist_temp=${j:1:3}

                if [ ${s} == 's' ]
                then
                        sim='stacked'
                        dim1='22.0'
                        dim2='7.0'
                        dim3='7.0'
                        pdb="../../pull_s.mdp"
                elif [ ${s} == 'p' ]
                then
                        sim='parallel'
                        dim1='16.0'
                        dim2='16.0'
                        dim3='7.0'
                        pdb="../../pull_p.mdp"
                elif [ ${s} == 'a' ]
                then
                        sim='antiparallel'
                        dim1='13.0'
                        dim2='16.0'
                        dim3='7.0'
                        pdb="../../pull_a.mdp"
                fi
                dist='s_'${dist_temp}
                pdb="pull.mdp"

                #bash Pull_sim_setup.bash ${file} ${sim} ${dist} ${pdb} ${dim1} ${dim2} ${dim3} 

                #bash Pull_sim_ndx_setup.bash ${file} ${sim} ${dist} ${pdb} ${dim1} ${dim2} ${dim3}     
                for l in {0..9}
                do
                        cd ~/AASims/inv_us/${file}/${sim}/${dist}/d_${l}
                        pwd
#bash ../../../get_distances.sh 
                        #bash ../../../setup_us_files.bash 
	   		cp ../../../npt_us.mdp .
			cp ../../../md_us.mdp .		
                done
#bash Pull_sim_equil.bash ${file} ${sim} ${dist} ${pdb} ${dim1} ${dim2} ${dim3} 
                #bash Pull_sim_pulling.bash ${file} ${sim} ${dist} ${pdb} ${dim1} ${dim2} ${dim3}       
        done
done

