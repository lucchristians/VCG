#!/bin/bash

#Variables
path=$1
id=$2
#solModel=$2
#temp=$3
#run=$4
#source /usr/local/gromacs-2021.1/bin/GMXRC
source /usr/local/gromacs-2022-plumed/bin/GMXRC

export GMX_FORCE_UPDATE_DEFAUT_GPU=true
export CUDA_VISIBLE_DEVICES=0,1,2
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
#export GMX_FORCE_UPDATE_DEFAULT_GPU=true
#export __GL_THREADED_OPTIMIZATIONS=0
#ulimit -c unlimited
#minimizing the energy of the syste
#if [ ! -f nvt.gro ]
#then
	gmx_mpi grompp -f nvt.mdp -c minimize.gro -r minimize.gro -p topol.top -o nvt.tpr -maxwarn 1
	gmx_mpi mdrun -v -deffnm nvt -nb gpu -bonded gpu -pme gpu -update cpu -pin auto -gpu_id ${id} 
#fi
#cd ..
#if [ ! -f npt.gro ] && [ -f nvt.gro ]
#then
	gmx_mpi grompp -f npt_us.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
	gmx_mpi mdrun -v -deffnm npt -nb gpu -bonded gpu -pme gpu -update cpu -pin auto -gpu_id ${id}
#fi
