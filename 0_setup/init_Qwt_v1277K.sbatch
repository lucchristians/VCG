#!/bin/bash

#SBATCH --job-name=i_Qwt_v1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100:4
#SBATCH --export=ALL
#SBATCH -p apak
#SBATCH --time=24:00:00
#SBATCH --mail-user=lucchristians@mines.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o gromacs-%j.out
#SBATCH -e gromacs-%j.err
#SBATCH -A 2104080412

#export GMX_GPU_DD_COMMS=true #only needed for multi-GPU
#export GMX_GPU_PME_PP_COMMS=true #only needed for multi-GPU

export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export CUDA_VISIBLE_DEVICES=0,1,2,3

##########################################################
# SUBMIT THIS JOB FROM HOME
# BUT PREPARE AND RUN THE JOB IN SCRATCH

# WE WILL MAKE 4 INSTANCES
cd $SCRATCH

#cp -r $HOME/newt/new/6xr8 ./new
SINGULARITY="singularity run --nv /sw/apps/singularity-images/gromacs-2021.sif"
DirA="Qwt_v1/tip3p/277K/run0"
DirB="Qwt_v1/tip3p/277K/run1"
DirC="Qwt_v1/tip3p/277K/run2"
DirD="Qwt_v1/tip3p/277K/run3"


cd $SCRATCH/new
for ((i=0; i<=3; i++ ))
do
	if [ $i -eq 0 ] 
	then
		cd $SCRATCH/new/${DirA}
	elif [ $i -eq 1 ]
	then
		cd $SCRATCH/new/${DirB}
	elif [ $i -eq 2 ]
	then
		cd $SCRATCH/new/${DirC}
	elif [ $i -eq 3 ]
	then
		cd $SCRATCH/new/${DirD}
	fi
	${SINGULARITY} gmx grompp -f nvt.mdp -o nvt.tpr -c minimize.gro -r minimize.gro -p topol.top -maxwarn 1

	let pinoff=$i*8
	
	#if [ ! -f nvt.cpt ]
	#then
    		${SINGULARITY} gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -pinoffset ${pinoff} -gpu_id ${i} &
	#else
    	#	${SINGULARITY} gmx mdrun -deffnm nvt -cpi nvt.cpt -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -pinoffset ${pinoff} -gpu_id ${i} &
	#fi
	cd $SCRATCH/new
done
wait	

for ((i=0; i<=3; i++ ))
do
	if [ $i -eq 0 ] 
	then
		cd $SCRATCH/new/${DirA}
	elif [ $i -eq 1 ]
	then
		cd $SCRATCH/new/${DirB}
	elif [ $i -eq 2 ]
	then
		cd $SCRATCH/new/${DirC}
	elif [ $i -eq 3 ]
	then
		cd $SCRATCH/new/${DirD}
	fi
	${SINGULARITY} gmx grompp -f npt.mdp -o npt.tpr -c nvt.gro -r nvt.gro -p topol.top -maxwarn 1

	let pinoff=$i*8
	
	#if [ ! -f nvt.cpt ]
	#then
    		${SINGULARITY} gmx mdrun -deffnm npt -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -pinoffset ${pinoff} -gpu_id ${i} &
	#else
    	#	${SINGULARITY} gmx mdrun -deffnm nvt -cpi nvt.cpt -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu -pme gpu -update gpu -pin on -pinoffset ${pinoff} -gpu_id ${i} &
	#fi
	cd $SCRATCH/new
done
wait	

##########################################################
