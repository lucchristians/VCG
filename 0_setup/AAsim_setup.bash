#!/bin/bash
#Variables
file=$1 #system name
pdb=$2 #protein structure file
dim1=$3 #nm 
dim2=$4 #nm
dim3=$5 #nm
path=$(pwd)

# Path to local gromacs installation
source /usr/local/gromacs-2021.1/bin/GMXRC

# Set up directories such that your system name is your file name
cd ${path}/${file}
# Convert to gro file type for appropriate model
echo 1 | gmx pdb2gmx -f ${pdb} -o  ${file}.gro -water tip3p -ignh  #ignore hydrogen atoms


# Determine simulation environment topology
echo 1 | gmx editconf -f ${file}.gro -o ${file}_bound.gro -bt triclinic -box ${dim1} ${dim2} ${dim3} -c

# Solvate using appropriate solvent model files
gmx solvate -cp ${file}_bound.gro -cs spc216.gro -o ${file}_Solv.gro -p topol.top

# Nuetralize the charge of the system (make sure to set up the ionize.mdp file to suit the simulation
gmx grompp -f ionize.mdp -c ${file}_Solv.gro -o ionize.tpr -p topol.top -maxwarn 1

# Ionize using sodium atoms for positive charge and chloride iones for negative charge
echo 13 | gmx genion -s ionize.tpr -o ${file}_IONIZE.gro -p topol.top -pname NA -nname CL  -neutral -conc 0.500

# Inital Equilibrations
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export __GL_THREADED_OPTIMIZATIONS=0
ulimit -c unlimited

# Minimizing the energy of the syste
gmx grompp -f minimize.mdp -c ${file}_IONIZE.gro -r ${file}_IONIZE.gro -p topol.top -o minimize.tpr -maxwarn 1
gmx mdrun -v -deffnm minimize -ntmpi 1 -ntomp 8
