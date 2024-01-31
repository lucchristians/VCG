#!/bin/bash               

#SBATCH -J ASSEM_VS
##SBATCH -o %j.out
##SBATCH -e %j.err
#SBATCH -p standard
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --mem-per-cpu=1800MB
#SBATCH -t 24:00:00
#SBATCH --mail-user=lucchristians@mines.edu
#SBATCH --mail-type=end,fail
#SBATCH -A bio220015

module load mvapich2/2.3.6
module load python
module load anaconda
module use $HOME/privatemodules
module load conda-env/mscg-py3.8.8
LMP=/anvil/projects/x-bio220015/lammps_bin/lmp_anvil

name=$1

#copy sim to scratch
for i in {1..4}
do
cp -r ${name} $SCRATCH/fin_assemb/${name}_${i}
done

#run replicas of the sim 
for i in {1..4}
do
        cd $SCRATCH/fin_assemb/${name}_${i}
        srun --exclusive -n 32 --mpi=pmi2 $LMP -i input_a -l log_${i}.txt -var SEED ${RANDOM} -var ITER 1 > /dev/null 2>&1 &
done
wait

