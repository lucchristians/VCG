#!/bin/bash

# Add virtual sites
#python lammps2pdb.py CA.lammpstrj 
python vs_setup.py CA.lammpstrj CA.pdb vs.lammpstrj #TODO: update this script

# Create virtual site containing pdb file
python lammps2pdb.py vs.lammpstrj
mv vs.pdb vs_ss.pdb
vmd -e alignment_CC.tcl -dispdev text

conda activate py3_7 #environment containing openMSCG
bash merge.bash

# Create henm files 
bash part0_create_sh.bash
bash part1_run_cghenm.bash

# TODO update this to include all processign
cp rcut_15.0/result.txt .

python create_vs_lt.py cg_mass.dat cg_charge.dat vs_frame.pdb result.txt vs_loc_stacked --outlt outfile.lt --gauss -1
~/progs/moltemplate_2021-5-17/moltemplate/scripts/moltemplate.sh outfile.lt -pdb vs_ss.pdb


