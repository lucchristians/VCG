echo "CG"
python rmsd_calc.py gauss/aligned_m.lammpstrj gauss/AA_m.pdb
echo "VCG"
python rmsd_calc.py vsite/aligned_m.lammpstrj vsite/vs_frame_m.pdb
