#!/bin/bash
#Variables
path=$1
id=$2
pdb='pull.mdp'
#solModel=$2
#temp=$3
#run=$4
#source /usr/local/gromacs-2021.1/bin/GMXRC
source /usr/local/gromacs-2022-plumed/bin/GMXRC
#if [ ! -z "${run}" ] 
#then
#	temps=("${temp}/${run}")
#elif [ ! -z ${temp} ]
#then
#	temps=(${temp})
#else
#	temps=(*K/)	
#fi
cd ${path}
#Inital Equilibrations
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export __GL_THREADED_OPTIMIZATIONS=0
ulimit -c unlimited
#minimizing the energy of the syste
gmx_mpi grompp -f ${pdb} -c npt.gro -r npt.gro -p topol.top -n sys.ndx -o pull.tpr -maxwarn 1
gmx_mpi mdrun -v -deffnm pull -plumed us_test.dat -ntomp 8 -nb gpu -bonded gpu -pme gpu -update gpu -pin auto -gpu_id ${id}
#cd ..
#gmx grompp -f ../../npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
#gmx mdrun -v -deffnm npt -nb gpu -bonded gpu -pme gpu -update cpu -pin on 
