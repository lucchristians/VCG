#!/bin/bash               

#SBATCH -J REM1
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
LAMMPS=/anvil/projects/x-bio220015/lammps_bin/lmp_anvil

#Directories
name1=rem_reps_novs_stacked
name2=rem_reps_novs_stacked_nointra
name3=rem_reps_novs_stacked_bc
name4=rem_reps_novs_stacked_nointra_bc

#copy directory to scratch x times
for i in {1..4}
do
	cp -r ${name1} $SCRATCH/${name1}_${i}
done

for i in {1..4}
do
	cp -r ${name2} $SCRATCH/${name2}_${i}
done

for i in {1..4}
do
	cp -r ${name3} $SCRATCH/${name3}_${i}
done

for i in {1..4}
do
	cp -r ${name4} $SCRATCH/${name4}_${i}
done

novx_iter=60
# run jobs
cd $SCRATCH/${name1}_1
bash rem_iterate.bash $LAMMPS -7.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name1}_2
bash rem_iterate.bash $LAMMPS -2.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name1}_3
bash rem_iterate.bash $LAMMPS -10.00 1.0 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name1}_4
bash rem_iterate.bash $LAMMPS -8.00 0.2 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/


cd $SCRATCH/${name2}_1
bash rem_iterate.bash $LAMMPS -7.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name2}_2
bash rem_iterate.bash $LAMMPS -2.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name2}_3
bash rem_iterate.bash $LAMMPS -10.00 1.0 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name2}_4
bash rem_iterate.bash $LAMMPS -8.00 0.2 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/

cd $SCRATCH/${name3}_1
bash rem_iterate.bash $LAMMPS -7.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name3}_2
bash rem_iterate.bash $LAMMPS -2.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name3}_3
bash rem_iterate.bash $LAMMPS -10.00 1.0 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name3}_4
bash rem_iterate.bash $LAMMPS -8.00 0.2 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/

cd $SCRATCH/${name4}_1
bash rem_iterate.bash $LAMMPS -7.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name4}_2
bash rem_iterate.bash $LAMMPS -2.00 0.6 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name4}_3
bash rem_iterate.bash $LAMMPS -10.00 1.0 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
cd $SCRATCH/${name4}_4
bash rem_iterate.bash $LAMMPS -8.00 0.2 8 ${novx_iter} & 
cd ~/anvil_rem_rep_stacked/
wait
