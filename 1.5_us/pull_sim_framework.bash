#!/bin/bash
#Variables
for i in "Qwt" 
do
	for j in "p0"
	do
		basedir="$(pwd)"
		cd ${basedir}
		file=${i}'_us'
		s=${j:0:1}
		dist_temp=${j:1:3}
		pdb="${i}_us.pdb"
		sim='parallel'
		dim1='12.0'
		dim2='17.0'
		dim3='9.0'
		dat="pull.mdp"
		
		dist='s_'${dist_temp}
		#create relivant pdbs
		python us_pdb_setup.py ${file} ${s} ${dist_temp}
		#set up directory structure
		bash pull_directory_setup.bash ${file} ${sim} ${dist} ${pdb} ${s} ${dim2} ${dim3} ${basedir}
		#set up simulation files
		bash Pull_sim_setup.bash ${basedir}/${file}/ ${pdb} ${file} ${dim1} ${dim2} ${dim3}	
		#create index file to define pull groups
		bash Pull_sim_ndx_setup.bash ${basedir}/${file}/	
		#setups up restraint and adjusts topol files based on common rule
		python pull_lock_setup.py ${basedir}/${file}/ 1 'A'
		

		# runs nvt and npt equilibrations
		bash Pull_sim_equil.bash ${file} 2
		# runs the pulling simulation
		bash Pull_sim_pulling.bash ${file} 2	
		
	done
done

