import numpy as np
import mdtraj as md
from mscg import *

#get center of mass from xyz input
def calc_com(xyz): #xyz [=] {atoms,3}
    return np.mean(xyz,axis=0)

#read in trajectory
traj=Trajectory('AA.lammpstrj', fmt='lammpstrj')
attr=Trajectory('AA_repl.lammpstrj', fmt='lammpstrj', mode='w')
ref_atoms=[27,34,41,48]

while(traj.read_frame()):
    atomtype=traj.t
    pos=traj.x
    attr.box=traj.box
    #TODO: redefine types to include the extra atoms
    newtypes=[]
    new_pos=[]
    max_type=atomtype[-1]
    chains=len(atomtype[max_type==atomtype])
    attr.natoms=traj.natoms+len(ref_atoms)*chains
    for j in range(chains):
        pent=int(j/5)
        newtypes.append(atomtype[j*max_type:(j+1)*max_type])
        new_pos.append(pos[j*max_type:(j+1)*max_type,:])
        for i in range(len(ref_atoms)):
            newtypes.append(np.array([max_type+i+1]))
            ref=np.copy(pos[atomtype==ref_atoms[i],:])
            com=np.copy(calc_com(pos[atomtype==ref_atoms[i],:][5*pent:5*(pent+1)]))
            new=(((com-ref[j])*0.45)+ref[j])
            
            new_pos.append(np.copy(new))
    newtypes=np.hstack(newtypes)
    new_pos=np.vstack(new_pos)
    attr.t=newtypes
    attr.x=new_pos
    attr.write_frame()
