#!/bin/bash
source /usr/local/gromacs-2021.3/bin/GMXRC

cd /mnt/data0/trrstore/barbar/tip3p/298K

# Defining selections for the top and bottom coiled coils
gmx make_ndx -f minimize.tpr -o sys.ndx < CC_index_input.txt

# CG sites are mapped to CA positions
echo 3 | gmx editconf -f minimize.gro -n sys.ndx -o CA.gro

for j in "run0" "run1" "run2" "run3"
do
	cd ${j}
	# Remove solvent in equilibration simulations
	echo 1 | gmx trjconv -f nvt.trr -n ../sys.ndx -skip 10 -o nvt.xtc
	echo 1 | gmx trjconv -f npt.trr -n ../sys.ndx -skip 10 -o npt.xtc
	
	gmx trjcat -f nvt.xtc npt.xtc step?.xtc step??.xtc -cat -o concatenated.xtc
	
	# This trajectory can be converted to a lammpstrj trajectory using vmd
	echo 3 | gmx trjconv -f concatenated.xtc -n sys.ndx -o concatenated_CA.xtc
	#
	echo 0 | gmx trjconv -f concatenated_CA.xtc -s CA.gro -o nojump_CA.xtc  -pbc nojump
	
	# aligns your molecule to the reference molecule (here, subset.pdb)
	echo 0 0 | gmx trjconv -f nojump_CA.xtc -s CA.gro -o aligned_CA.xtc -fit rot+trans -force yes -vel yes
	echo 19 19 | gmx trjconv -f aligned_CA.xtc -n sys.ndx -s CA.gro -o aligned_top.xtc -fit rot+trans -force yes -vel yes
	echo 20 | gmx trjconv -f aligned_CA.xtc -n sys.ndx -s CA.gro -o aligned_bot.xtc -force yes -vel yes
	echo 19 19 | gmx trjconv -f aligned_bot.xtc -n sys.ndx -s CA.gro -o aligned_bot_t.xtc -fit rot+trans -force yes -vel yes
	echo 19 | gmx rms -f aligned_top.xtc -s CA_1.gro -n sys.ndx -o rmsd_top.xvg -fit none
	echo 19 | gmx rms -f aligned_bot_t.xtc -s CA_1.gro -n sys.ndx -o rmsd_bot.xvg -fit none
	cd ..
done
