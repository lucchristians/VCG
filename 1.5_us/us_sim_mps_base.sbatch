#!/bin/bash

#SBATCH --job-name={name}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4 
#SBATCH --gres=gpu:a100:1
#SBATCH --export=ALL
#SBATCH -p apak
#SBATCH --time=24:00:00
#SBATCH --mail-user=lucchristians@mines.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o gromacs-%j.out
#SBATCH -e gromacs-%j.err
#SBATCH -A 2104080412
#SBATCH -C mps
#SBATCH --mem=4693M

#export GMX_GPU_DD_COMMS=true #only needed for multi-GPU
#export GMX_GPU_PME_PP_COMMS=true #only needed for multi-GPU
#export GMX_FORCE_UPDATE_DEFAULT_GPU=true
#export CUDA_VISIBLE_DEVICES=0,1,2,3

################################################

nvidia-cuda-mps-control -d
cd $SCRATCH
SINGULARITY="singularity run --nv /sw/apps/singularity-images/gromacs_2021-dev-plumed.sif"
DirA={path}
DirB={path}_5
cd $SCRATCH/${DirA}
DIR=${DirA}
${SINGULARITY} gmx grompp -f ../md_us.mdp -o step2.tpr -c step1.gro -r step1.gro -p topol.top -maxwarn 1
${SINGULARITY} gmx mdrun -deffnm step2 -plumed us.dat -ntomp 4 -nb gpu -bonded gpu -pme gpu -update cpu -pin off &
cd $SCRATCH/${DirB}
DIR=${DirB}
${SINGULARITY} gmx grompp -f ../md_us.mdp -o step2.tpr -c step1.gro -r step1.gro -p topol.top -maxwarn 1
${SINGULARITY} gmx mdrun -deffnm step2 -plumed us.dat -ntomp 4 -nb gpu -bonded gpu -pme gpu -update cpu -pin off &
wait

################################################
