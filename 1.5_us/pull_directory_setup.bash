#!/bin/bash
#Variables
file=$1
sim=$2
dist=$3
pdb=$4
dim1=$5
dim2=$6
dim3=$7
path=$8

if [ ! -d ${file} ]
then
	mkdir ${file}
fi
cd ${file}
if [ ! -d ${sim} ]
then
	mkdir ${sim}
fi
cd ${sim}
if [ ! -d ${dist} ]
then
	mkdir ${dist}
fi
cd ${dist}
#TODO: move mdp and ff directorys 
#pwd
cp ../../../*.mdp .
mv "pull_${dim1}.mdp" "pull.mdp"
#pullInds=( $(python ../../../define_pull_atoms.py) )
#sed -i -z "s/${pullInds[0]}/{atom1}/g" pull.mdp
#sed -i -z "s/${pullInds[1]}/{atom2}/g" pull.mdp
cp -r ../../../*.ff .

