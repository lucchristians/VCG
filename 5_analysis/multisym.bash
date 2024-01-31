echo "CG"
python sym_calc.py gauss/aligned_m.lammpstrj gauss/AA_m.pdb sym5_CG.txt
echo "VCG"
python sym_calc.py vsite/aligned_m.lammpstrj vsite/vs_frame_m.pdb sym5_VCG.txt
