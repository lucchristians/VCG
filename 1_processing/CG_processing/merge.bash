rm vs_m.lammpstrj
for (( i=1; i<5; i++ ))
do
	cat CA_monomers${i}.lammpstrj >> vs_m.lammpstrj
	#rm CA_monomers${i}.lammpstrj
done

