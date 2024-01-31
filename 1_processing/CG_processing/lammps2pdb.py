import numpy as np
import sys as sys

#inputs are either:
#   lammpstrj_file
#   lammpstrj_file #max_#_of_chains
print(len(sys.argv))
if(len(sys.argv) != 2 and len(sys.argv) != 3):
    raise Exception("incorrect number of arguments (takes 1 or 2 arguments)")

traj_in = sys.argv[1]

if( len(sys.argv) == 3):
    max_chains = sys.argv[2]
else:
    max_chains = 9999
print(max_chains)
traj = open(traj_in,'r')
pdb = open(traj_in[:traj_in.find('.')]+".pdb",'w')
for i in range(9):
    traj.readline()
line= traj.readline()
params = line.split()
chain = 0
chain_prev = 0 
index = 0
residue = 0
residue_prev = 0
pos=[]
while(not (params[0] == "ITEM:")):
    index = params[0]
    atom_name = "CG"
    residue = params[1]
    resname = "C%-2d"%(int(residue))
    if(int(residue_prev) > int(residue)):
        chain+=1
    if(int(chain)>=int(max_chains)):
        print("limit reached")
        break
    pos1=float(params[2])
    pos2=float(params[3])
    pos3=float(params[4])
    input_line="%4s  %5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%("ATOM",int(index),atom_name,resname,chain,int(residue),pos1,pos2,pos3,0.00,0.00)
    pdb.write(input_line)
    residue_prev=residue
    line = traj.readline()
    params = line.split()
    if(len(params)==0):
        break

pdb.write("END\n")
pdb.close()
traj.close()
