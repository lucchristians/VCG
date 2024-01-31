#!/bin/bash
#Variables
path=$1
pdb=$2
file=$3
dim1=$4
dim2=$5
dim3=$6
source /usr/local/gromacs-2021.1/bin/GMXRC

cd ${path}
#convert to gro file type for appropriate model

#for i in "${temps[@]}"
#do
#	cd ${i}
	echo 1 | gmx pdb2gmx -f ${pdb} -o  ${file}.gro -water tip3p -ignh  #ignore hydrogen atoms

	pullInds=( $(python ../define_pull_atoms.py ${file}.gro) )
	#sed -i -z "s/{atom1}/${pullInds[0]}/g" pull.mdp
	#sed -i -z "s/{atom2}/${pullInds[1]}/g" pull.mdp
	#sed -i -z "s/{atom3}/${pullInds[2]}/g" pull.mdp
	#sed -i -z "s/{atom4}/${pullInds[3]}/g" pull.mdp

	#determine simulation environment topology
	#echo 1 | gmx editconf -f ${file}.gro -o ${file}_bound.gro -bt triclinic -box 15.00 11.00 6.50 -c -princ
	#echo 1 | gmx editconf -f ${file}.gro -o ${file}_bound.gro -bt triclinic -box ${dim1} ${dim2} ${dim3} -c -princ
	echo 1 | gmx editconf -f ${file}.gro -o ${file}_bound.gro -bt triclinic -box ${dim1} ${dim2} ${dim3} -c
	#echo 1 | gmx editconf -f ${file}.gro -o ${file}_bound.gro -bt triclinic -d 2.5 -c -princ
	
	#solvate using appropriate solvent model files
	gmx solvate -cp ${file}_bound.gro -cs spc216.gro -o ${file}_Solv.gro -p topol.top
	
	#nuetralize the charge of the system (make sure to set up the ionize.mdp file to suit the simulation
	gmx grompp -f ionize.mdp -c ${file}_Solv.gro -o ionize.tpr -p topol.top -maxwarn 1
	# ionize using sodium atoms for positive charge and chloride iones for negative charge
	echo 13 | gmx genion -s ionize.tpr -o ${file}_IONIZE.gro -p topol.top -pname NA -nname CL  -neutral -conc 0.500
	
	#Inital Equilibrations
	export GMX_FORCE_UPDATE_DEFAULT_GPU=true
	export __GL_THREADED_OPTIMIZATIONS=0
	ulimit -c unlimited
	#minimizing the energy of the syste
	gmx grompp -f minimize.mdp -c ${file}_IONIZE.gro -r ${file}_IONIZE.gro -p topol.top -o minimize.tpr -maxwarn 1
	gmx mdrun -v -deffnm minimize
	#cd ..
