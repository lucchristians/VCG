#!/bin/bash
#Variables
path=$1
#solModel=$2
#temp=$3
#run=$4
source /usr/local/gromacs-2022.5/bin/GMXRC

cd ${path}
#Inital Equilibrations
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export __GL_THREADED_OPTIMIZATIONS=0
ulimit -c unlimited
#minimizing the energy of the syste
gmx make_ndx -f minimize.tpr -o sys.ndx < ../sys_input.txt 
