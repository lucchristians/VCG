rm henm_params_ref.sh

#awk '
#BEGIN{
#printf "--mass ";
#}{
#if(NR>1) {
#printf "%f ",$2;
#}
#}
#END{printf "\n"}' < cg_mass.dat >> henm_params_ref.sh

echo "--traj=../imflex_a.lammpstrj" >> henm_params_ref.sh
echo "-v 9" >> henm_params_ref.sh
echo "--cut (CUT)" >> henm_params_ref.sh
echo "--ktol 0.01" >> henm_params_ref.sh
echo "--sdftol 0.001" >> henm_params_ref.sh
echo "--save result" >> henm_params_ref.sh
