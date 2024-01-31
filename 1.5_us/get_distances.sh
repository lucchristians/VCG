#!/bin/bash

#################################################
# get_distances.sh
#
#   Script iteratively calls gmx distance and
#   assembles a series of COM distance values
#   indexed by frame number, for use in
#   preparing umbrella sampling windows.
#
# Written by: Justin A. Lemkul, Ph.D.
#    Contact: jalemkul@vt.edu
#
#################################################
source /usr/local/gromacs-2022.5/bin/GMXRC
path=$1
direction=$2 #2-x 3-y 4-z
coord="{print\$${direction}}"
cd ${path}
echo 0 | gmx trjconv -s pull.tpr -f pull.trr -o conf.gro -sep

# compute distances
for (( i=0; i<51; i++ ))
do
    gmx distance -s pull.tpr -f conf${i}.gro -n sys.ndx -select 'com of group "REST" plus com of group "PULLT"' -oall dist${i}.xvg -oxyz distxyz${i}.xvg 
done

# compile summary

rm summary_distances.dat
rm summary_distancesxyz.dat
touch summary_distances.dat
touch summary_distancesxyz.dat

for (( i=0; i<51; i++ ))
do
    d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
    echo "${i} ${d}" >> summary_distances.dat
    rm dist${i}.xvg
    d=`tail -n 1 distxyz${i}.xvg | awk ${coord}`
    echo "${i} ${d}" >> summary_distancesxyz.dat
    rm distxyz${i}.xvg
done

exit;
