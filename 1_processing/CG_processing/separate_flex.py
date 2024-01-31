import numpy as np
from mscg import *
import sys

def fix_pbc(xyz,unitcell):
    xyz_=np.copy(xyz)
    for m in range(3):
        #determine how out of the box
        pbc_shift = np.floor(xyz[:,:,m]/unitcell[m])
        xyz_[:,:,m] = xyz_[:,:,m] - pbc_shift*unitcell[m]
    return np.copy(xyz_)

def shift_com(xyz,pos,unitcell):
    xyz_=np.copy(xyz)
    com=np.sum(xyz_,axis=1)/np.shape(xyz_)[1]
    shift=-com+pos
    xyz_shifted = xyz_+np.expand_dims(shift,axis=1)
    return fix_pbc(xyz_shifted,unitcell)
			
#TODO: write lammps frames here
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

def write_lammps_frame(fopen,num,xyz,unitcell,types): #assumes each atom is of a unique type and of one chain
    n_atoms=np.shape(xyz)[0]
    write_lammps_header(fopen,num,n_atoms,unitcell)
    for i in range(n_atoms):
        print_atom(fopen,i,types[i],xyz[i,:])


traj_name = sys.argv[1] #format lammpstrj
# Load in, reads & process trajectory
print("load traj")
traj=Trajectory(traj_name, fmt='lammpstrj')
pos_array = []
while(traj.read_frame()):
    atomtype=traj.t
    unitcell=traj.box
    pos=traj.x
    pos_array.append(np.copy(pos))
n_frames=len(pos_array)
pos_array=np.array(pos_array)
print("centering atoms")
#split trajectory
flex_region=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,54])
imflex_region=np.arange(17,59,dtype=int)
flex_=pos_array[:,flex_region,:]
flex=shift_com(flex_,unitcell/3,unitcell)
flex=shift_com(flex,unitcell/2,unitcell)
flex=shift_com(flex,unitcell/2,unitcell)
del flex_
flex_t=atomtype[flex_region]

imflex_=pos_array[:,imflex_region,:]
imflex=shift_com(imflex_,unitcell/3,unitcell)
imflex=shift_com(imflex,unitcell/2,unitcell)
imflex=shift_com(imflex,unitcell/2,unitcell)
del imflex_
print("creating files")
imflex_t=atomtype[imflex_region]
f_flex=open("flex.lammpstrj",'w')
f_imflex=open("imflex.lammpstrj",'w')
for i in range(len(pos_array)):
    write_lammps_frame(f_flex,i,flex[i],unitcell,flex_t)
    write_lammps_frame(f_imflex,i,imflex[i],unitcell,imflex_t)

f_flex.close()
f_imflex.close()
