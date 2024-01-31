#!/bin/bash
path=$1
cd ${path}
files=( $(python ../setup_us_files.py) )

for i in {0..29}
do
	file="${files[${i}]}"
	mkdir "d_${i}"
	cd "d_${i}"

	cp ../conf${file}.gro .
	mv conf${file}.gro d_${i}.gro
	cp -r ../charmm36-feb2021.ff .
	cp ../topol* .
	cp ../posre* .
	cp ../*.mdp .

	cd ..

done
