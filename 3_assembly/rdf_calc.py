import numpy as np
import networkx as nx
import mdtraj as md
import sys as sys
import math as math
from mscg import *
# USAGE:
# analyze_assembly_and_extract_clusters.py trajbasename startindex endindex clusterMolID dtau_per_frame 
# this version assumes trajectory is LAMMPSTRJ format custom (x, y, z)

##################### USER SPECIFIED INFO ##############################
path=sys.argv[1]
if(path[-1]!='/'):
    path+='/'
trajname=sys.argv[2] # The name of the trajectory (base+index.lammpstrj)
#pdbname=sys.argv[3]
step_id=int(sys.argv[3])
##################### END USER SPECIFIED INFO ##############################

#TODO: modify parameters

#TODO: modify to reflect desired sites
#pi_bonders={'HIS':['SC2','SC3']}
#TODO: verify this indexing
pairs=[[22,23],[29,30],[36,37],[43,44],[50,51]]

#returns array ( nframes, natoms, (aid,type,x,y,z) )
def load_atom_traj(trajname) : 
    i=[]
    t=[]
    pos_array = []
    #print(i)
    traj=Trajectory(trajname, fmt='lammpstrj')
    num=-1
    while(traj.read_frame()):
        num+=1
        #if(num%10==0):
        atomtype=traj.t
        box=traj.box
        i.append(np.arange(0,len(atomtype),1).astype(int).reshape(-1,1))
        t.append(atomtype.reshape(-1,1))
        pos_array.append(np.copy(traj.x))
    nframes=len(pos_array)
    natoms=len(atomtype)
    i=np.array(i)
    t=np.array(t)
    pos_array=np.array(pos_array)
    #it=np.append(i,t,axis=2)
    #traj_array=np.append(it,pos_array,axis=2)
    
    return i, t, pos_array, box, nframes, natoms

def load_atom_traj_trr(trajname,topo) : 
    traj = md.load(trajname,top=topo)
    nframes=traj.n_frames
    box=traj.unitcell_lengths
    natoms=traj.n_atoms
    return traj, box, nframes, natoms

def fast_dist(traj,box,ind_i,ind_j):
    dist=[]
    new_box=box.reshape([1,1,-1])
    len_j=np.arange(len(ind_j))
    for i in range(len(ind_i)):
        ii=ind_i[i]
        jj=ind_j[len_j!=i]
        xyz_i = traj[:,[ii],:]
        xyz_j = traj[:,jj,:]
        dr = xyz_i - xyz_j
        dr = dr - new_box * np.rint( dr / new_box )
        dist.append(np.linalg.norm(dr,axis=-1))
    return np.hstack(dist)

def calc_dist(r1, r2, dims) :
    dr = r2 - r1
    dr = dr - dims * np.rint( dr / dims )
    dist = np.linalg.norm(dr,axis=-1)
    return dist

##################### COLLECT ALL TRAJECTORY DATA FIRST ######################################
nframes = 0
i,t,traj, dims, nframes, natoms = load_atom_traj(path+"%s.lammpstrj" % trajname) 
ncg=int(natoms/484)
#print(ncg)
#traj, dims, nframes, natoms = load_atom_traj_trr(path+trajname,path+pdbname) 
#print(traj)
##################### ANALYZE ALL TRAJECTORY DATA ######################################
count=1
for pair in pairs:
    ind_i_list=[]
    ind_j_list=[]
    for i in range(484):
        ind_i_list.append(i*ncg+pair[0])
        ind_j_list.append(i*ncg+pair[1])
    ind_i_list=np.array(ind_i_list)
    ind_j_list=np.array(ind_j_list)
    dist=fast_dist(traj,dims,ind_i_list,ind_j_list)
    dist=dist.flatten()
    dist=dist[dist<20.0]
    np.save(f'{path}pair{count}_step{step_id}.npy',dist)
    count+=1
