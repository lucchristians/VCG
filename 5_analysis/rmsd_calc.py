import numpy as np
import mdtraj as md
import sys

target_traj=sys.argv[1]
target_top=sys.argv[2]

targ=md.load(target_traj,top=target_top)
refe=md.load(target_top)

selection=targ.top.select('residue 18 to 48')
print(len(selection))
rmsd=md.rmsd(targ,refe,atom_indices=selection)
rmsd=rmsd[rmsd<4.0] # removes outliers
r_m=np.mean(rmsd)
r_s=np.std(rmsd)
print("RMSD: ",r_m,"+/-",r_s,"nm")

