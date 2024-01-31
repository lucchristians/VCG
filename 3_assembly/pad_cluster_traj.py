import numpy as np
import sys as sys
import math as math

# USAGE:
# pad_cluster_traj.py trajname 
# assumes trajectory is CUSTOM LAMMPSTRJ format (id, type, x, y, z, vx, vy, vz)

# Will find maximum number of atoms over all frames  
# Then pad each frame accordingly and correct atom indices (vx = 1) 
# (fake particles will have vx = 0 and x,y,z = 0)

##################### USER SPECIFIED INFO ##############################

trajname=sys.argv[1] # The name of the trajectory
ncg=int(sys.argv[2]) # number of sites per monomer (assumes only one type of monomer)

##################### END USER SPECIFIED INFO ##############################

# write lammps stuff
def write_lammps_header(tfile, time, natoms, dims, fields) :
    tfile.write("ITEM: TIMESTEP\n")
    tfile.write("%d\n" % time)
    tfile.write("ITEM: NUMBER OF ATOMS\n")
    tfile.write("%d\n" % natoms)
    tfile.write("ITEM: BOX BOUNDS pp pp pp\n")
    tfile.write("0.0 %f\n" % dims[0])
    tfile.write("0.0 %f\n" % dims[1])
    tfile.write("0.0 %f\n" % dims[2])
    tfile.write("ITEM: ATOMS ")
    for f in fields :
        tfile.write("%s " % f)
    tfile.write("\n")

def write_lammps_atom(tfile, index, type, x, y, z, vx, vy, vz) :
    tfile.write("%d %d %f %f %f %f %f %f\n" % (index, type, x, y, z, vx, vy, vz))


##################### COLLECT NMOL STATS OVER ALL TRAJECTORY DATA ######################################

FTRAJ=open(trajname,"r")
all_lines = FTRAJ.readlines()
ncg_max = 0

ncg_list = [] #number of cg sites per frame

dims = []

start_line = 0
fcount = 0
while(start_line < len(all_lines)) :
    
    natoms = int(all_lines[start_line+3].split()[0])
    print(natoms)
    [xlo, xhi] = [float(i) for i in all_lines[start_line+5].split()] 
    [ylo, yhi] = [float(i) for i in all_lines[start_line+6].split()] 
    [zlo, zhi] = [float(i) for i in all_lines[start_line+7].split()] 
    XLEN = xhi - xlo
    YLEN = yhi - ylo
    ZLEN = zhi - zlo
    dims = [XLEN, YLEN, ZLEN]

    start_line += 9

    print(f"Current frame has {natoms:d} atoms")
    ncg_list.append(natoms)
    if(natoms > ncg_max) :
        ncg_max = natoms
        
    start_line += natoms
    fcount += 1

nframes = len(ncg_list)
FTRAJ.close()

##################### PROCESS AND PAD ALL TRAJECTORY DATA ######################################
outfile = open("padded.lammpstrj","w")

FTRAJ=open(trajname,"r")
all_lines = FTRAJ.readlines()
start_line = 0
maxatoms = ncg_max

for i in range(nframes) :
    natoms = int(all_lines[start_line+3].split()[0])
    time = int(all_lines[start_line+1].split()[0])
    write_lammps_header(outfile, time, maxatoms, dims, ["id", "type", "x", "y", "z", "vx", "vy", "vz"])

    start_line += 9
    serial = 1

    # orig atoms
    for a in range(ncg_list[i]) :
        line = all_lines[start_line+a].split()
        aid = int(line[0])
        atype = int(line[1])
        x = float(line[2])
        y = float(line[3])
        z = float(line[4])
        fx = float(line[5])
        fy = float(line[6])
        fz = float(line[7])
        
        vx = 1 # real, 0 = fake
        vy = 0 
        vz = 0

        write_lammps_atom(outfile, serial, atype, x, y, z, vx, vy, vz)
        serial += 1

    # write pad atoms
    for m in range(ncg_max - ncg_list[i]) :
        atype = serial%ncg
        if(atype == 0) :
            atype = ncg
        write_lammps_atom(outfile, serial, atype, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        serial += 1
    
    start_line += ncg_list[i]


outfile.close()
FTRAJ.close()
