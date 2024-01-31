#rcut - cutoff distances in Angstroms
rcut=( 15.0 )
for r in "${rcut[@]}" 
do
    mkdir rcut_${r}
    cd rcut_${r}
    cp ../henm_params_ref.sh henm_params.sh
    sed -i "s/(CUT)/${r}/g" henm_params.sh
    cghenm @henm_params.sh
    cd ../
done


