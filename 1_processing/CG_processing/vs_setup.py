#!/usr/env/bin python
import mdtraj as md
import numpy as np
import sys as sys

#define what real site (based on mapped sites) each virual site will be based off of
vs_location=[24,31,38,45,52]
vs_location=np.array(vs_location)

#this script will assume the trajectory your using is either a lammpstrj with cg sites defined or you are using a CA-only trr trajectory
if len(sys.argv) != 4 :
    print("Usage: vs_setup.py traj_in top_in traj_out")
    sys.exit()

traj_in = sys.argv[1] # either trr or lammpstrj
traj_type = traj_in[traj_in.find('.')+1:]
top_in = sys.argv[2]
top_type = top_in[top_in.find('.')+1:]
    
traj_out = sys.argv[3]
    
###
# Defining functions
###
def write_lammps_header(fopen, time, natoms, unitcell) :
    fopen.write("ITEM: TIMESTEP\n")
    fopen.write("%d\n" % time)
    fopen.write("ITEM: NUMBER OF ATOMS\n")
    fopen.write("%d\n" % natoms)
    fopen.write("ITEM: BOX BOUNDS pp pp pp\n")
    fopen.write("%f %f\n" % (0.0, unitcell[0]))
    fopen.write("%f %f\n" % (0.0, unitcell[1]))
    fopen.write("%f %f\n" % (0.0, unitcell[2]))
    fopen.write("ITEM: ATOMS id type x y z\n")
    
def print_atom(fopen,site_ind,site_type,pos):
    fopen.write("%d %d %f %f %f\n"%(site_ind,site_type,pos[0],pos[1],pos[2]))

def get_cgsite_pos(selection,traj):
    pos = md.compute_center_of_mass(traj.atom_slice(traj.top.select(selection)))
    return pos


###
# script begins here
###

print("reading file...",end='\r',flush=True)
trj = md.load(traj_in, top=top_in)
top = trj.topology

f = open(traj_out, 'w')

#frame and box information
uc_lengths = trj.unitcell_lengths


print("mapping system...",end='\r')
n_cgsites = trj.n_atoms 

#defining virtual sites for one chain.

#determine residues for each cg site "based on map.txt"
#applys to lammpstrj
n_chains = top.n_chains
n_residues = top.chain(0).n_atoms

vs_sites=len(vs_location)*n_chains
sites=n_cgsites+vs_sites
sites_per_chain=sites/n_chains

cgsite_pos = [] #TODO: only need to do this if trrs are used as input
#lammpstrj case
ind=[]
for i in range(n_chains):
    ind.append(np.arange(i*n_residues,(i+1)*n_residues))
    ind.append((vs_location+(i+1)*n_residues)%(n_chains*n_residues))
ind=np.hstack(ind)
print(np.shape(ind[0]))
print(type(ind[0]))
for i in range(trj.n_frames):
    if(i%10==0):
        print("processing frame %7d out of %7d     "%(i,trj.n_frames),end='\r')
    write_lammps_header(f,i,len(ind),uc_lengths[i]*10.0)
    #print each site for a given frame
    for j in range(len(ind)):
        # lammpstrj case
        print_atom(f,j+1,(j%sites_per_chain)+1,trj.xyz[i,ind[j],:]*10.0)
            
         
print("done                                              ")


