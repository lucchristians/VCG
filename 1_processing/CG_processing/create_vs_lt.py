#!/usr/bin/env python

# This script will take the output files from EDCG and hENM:
# 
# 1) cg_mass.dat
# 2) cg_charge.dat
# 3) results.txt
# 4) vslocation.txt
#
# and the PDB file of the CG reference:
# 1) reference.pdb 
#
# to create:
# 1) A PDB file of the CG-mapped molecule
# 2) An LT file of the CG-mapped molecule (only excl vol and bonds)

import argparse
import sys
import os
import math
import numpy as np
from Bio.PDB import PDBParser, PDBIO

###
# Before use: define indexes (base 0) of vs_sites and the real sites they interact with. 
#             This can be done for vsites that interact with other vsites if defined twice
#             (i.e. '32 86' & '86 32')
###


parser = argparse.ArgumentParser(description='Reads the mapping and PDB files to create a mapped PDB and LT.')
#parser.add_argument('mapfile',metavar='MapFile', help='The name of the CG mapping file')
parser.add_argument('massfile',metavar='MassFile', help='The name of the CG mass file')
parser.add_argument('chargefile',metavar='ChargeFile', help='The name of the CG charge file')
parser.add_argument('pdbfile',metavar="PDBFile", help='The name of the reference CG PDB file')
parser.add_argument('henmfilei',metavar="HENMFile1", help='The name of the nvs hENM output file')
parser.add_argument('vslocation',metavar='VSLocation',help='maps each virtual site to an associated real site')
parser.add_argument('--outlt',metavar='LTOut', default='outfile.lt', help='The name of the output mapped LT file')
parser.add_argument('--gauss',metavar='gauss', default='3.0', help='the magnitude of the gauss pair style', type=str)
args = parser.parse_args()

#mapping = np.genfromtxt(args.mapfile, dtype=np.int)
masses = np.genfromtxt(args.massfile) #TODO: modify so it takes in mass and charge based on non virtual site model 
charges = np.genfromtxt(args.chargefile)
if(len(masses) != len(charges)) :
    print("Mismatch in the number of CG sites! Please correct input files\n")
    exit(1)

vs_index = np.genfromtxt(args.vslocation, dtype=int)
pdb = open(args.pdbfile,'r')
if(len(vs_index)>0):
    virt_sites_index=vs_index[:,1]
    real_sites_index=vs_index[:,0]
else:
    virt_sites_index=[]
    real_sites_index=[]

#obtaining neccessary information from the pdb file
xyz = [] #TODO: section is commented out in REM script; monomer number and chain number are explicitly defined
n_chains = 0
monomer_length=0
curr_id = -1
ignore_once=True
ignore_chain_id_difference=False
monomer_length_calculated=False
while(True): #iterate until end statement is met
    pdb_line = pdb.readline()
    #check if new chain has been reached
    if(pdb_line[0:3] == "TER"):
        n_chains += 1
        monomer_length_calculated=True
        ignore_chain_id_difference=True
    if(pdb_line[0:3] == "END"):
        break
    # assumes pdb has position data in angstroms TODO: correctly use nm in future
    xyz.append([float(pdb_line[30:38]),float(pdb_line[38:46]),float(pdb_line[46:54])]) #positions
    prev_id = curr_id
    curr_id = pdb_line[21] # chain id 
    #determine if a new chain is reached
    if(curr_id != prev_id and not ignore_chain_id_difference):
        n_chains += 1
        monomer_length_calculated=True 
    else:
        ignore_chain_id_difference=False
    if(not monomer_length_calculated):
        monomer_length += 1
    elif(ignore_once):
        monomer_length += 1
        monomer_length_calculated=False
        ignore_once=False 

print(n_chains,monomer_length)
#defining vs_sites mass and charge on real properties (charges are assumed zero
new_masses = []
new_charges = []
real_index=0
for i in range(monomer_length):
    vs_present=False
    for j in range(len(virt_sites_index)):
        if(virt_sites_index[j]==i):
            ind=j # to use the real site equivalent to define later
            
            #check for special case were both are virtual sites
            two_vs_present=False
            for k in range(len(real_sites_index)):
                if(real_sites_index[k] == virt_sites_index[j]):
                    two_vs_present=True
                    break
            vs_present=True
            break 
    if(vs_present):
        if(two_vs_present):
            new_masses.append([i+1,100.0])
            new_charges.append([i+1,0.0])
        else:
            new_masses.append([i+1,masses[real_sites_index[ind],1]])
            new_charges.append([i+1,0.0])
    else:
        new_masses.append([i+1,masses[real_index,1]])
        new_charges.append([i+1,charges[real_index,1]])
        real_index += 1    
charges = np.array(new_charges)
masses = np.array(new_masses)
pdb.close()

###
# seting up lt file below
###

#create CG-mapped LT file
LTfile = open(args.outlt, "w")
LTfile.write("CG {\n")

#atom block
LTfile.write(" write(\"Data Atoms\") {\n")
for i in range(monomer_length) : #TODO: change so mapping is not needed
    LTfile.write("   $atom:%d  $mol:1  @atom:CG%d  %f  %f  %f  %f \n" % (i+1, i+1, charges[i,1], xyz[i][0], xyz[i][1], xyz[i][2])) #TODO: positionwas were initially defined at origin in REM script
LTfile.write(" }\n")

#bond block
all_bonds = np.genfromtxt(args.henmfilei, skip_header=1)



#print bond information
LTfile.write(" write(\"Data Bonds\") {\n")
for i, bondpair in enumerate(all_bonds) :
    ii = bondpair[0]
    jj = bondpair[1]
    LTfile.write("   $bond:%d  @bond:%d_%d  $atom:%d  $atom:%d\n" % (i+1, ii, jj, ii, jj))
LTfile.write(" }\n")

#mass block
LTfile.write(" write_once(\"Data Masses\") {\n")
for i, m in enumerate(masses) :
    ii = m[0]
    mass = m[1]
    LTfile.write("   @atom:CG%d  %f\n" % (ii, mass))
LTfile.write(" }\n")

#Boundaries block
LTfile.write("write_once(\"Data Boundary\"){\n") #TODO: box sizes are different (not most important in implicit solvent
LTfile.write("00.0000 1200.0000 xlo xhi\n")
LTfile.write("00.0000 1200.0000 ylo yhi\n")
LTfile.write("00.0000 1200.0000 zlo zhi\n")
LTfile.write("}\n")

#settings block (pair/bond coeffs)
LTfile.write(" write_once(\"In Settings\") {\n")
#pair coeffs
for i in range(len(virt_sites_index)):
    if(float(args.gauss)<0.0):
    	LTfile.write("   pair_coeff  @atom:CG%d  @atom:CG%d  gauss  A%d B%d\n" % (virt_sites_index[i]+1,real_sites_index[i]+1,i+1,i+1))
    else:
    	LTfile.write("   pair_coeff  @atom:CG%d  @atom:CG%d  gauss  %5.2f  0.2\n" % (virt_sites_index[i]+1,real_sites_index[i]+1,float(args.gauss))) #TODO: based on constant placeholders 'An Bn' not explicit constants in REM script
    

for i in range(monomer_length) :
    for j in range(monomer_length) :
        ii = i+1
        jj = j+1
        vs_present=False #TODO: change to have exclusive 2 vs mode
        for k in range(len(virt_sites_index)):
            if(i==virt_sites_index[k] or j==virt_sites_index[k]): #TODO: chekc if a interaction is already defined
                vs_present=True
                break
        
        # avoid inverted interactions (i.e. 1-2 vs 2-1)
        if(i <= j) :
            if(vs_present):
                #exclude any virtual site soft exclusions
                if(i == j):
                    LTfile.write("   pair_coeff  @atom:CG%d  @atom:CG%d  soft  20.0000  4.5\n" % (ii,jj))
                else:
                    LTfile.write("   pair_coeff  @atom:CG%d  @atom:CG%d  soft  0.0000   4.5\n" % (ii,jj))
                continue
            LTfile.write("   pair_coeff  @atom:CG%d  @atom:CG%d  soft  20.0000   4.5\n" % (ii,jj))
            if(np.abs(charges[i,1])>0.0001 and np.abs(charges[j,1])>0.0001):
                LTfile.write("   pair_coeff  @atom:CG%d  @atom:CG%d  coul/debye\n" % (ii,jj))
LTfile.write("\n")
#bond coeffs
count_vs_pairs=1
k_const_ind=np.zeros(len(vs_index))
for i, bondpair in enumerate(all_bonds) :
    ii = bondpair[0]
    jj = bondpair[1]
    r0 = bondpair[2]
    k = bondpair[3] #TODO: added is vs part to here in REM script 
    # factor of 0.5 is because LAMMPS doesn't use 1/2*k*dr^2
    is_vs_pair=False 
    for j in range(len(vs_index)):
        if((ii == vs_index[j,0]+1 and jj == vs_index[j,1]+1) or (ii == vs_index[j,1]+1 and jj == vs_index[j,0]+1)):
            is_vs_pair=True
            break
    if(is_vs_pair and float(args.gauss) < 0.0):
        LTfile.write("   bond_coeff  @bond:%d_%d  K%d   %f\n" % (ii, jj, j+1, r0))
        k_const_ind[j]=i+1
        #k_const_ind.write("%d\n"%(i+1))
        count_vs_pairs+=1 #TODO: make work for 2 vs 
    else:
        LTfile.write("   bond_coeff  @bond:%d_%d  %f   %f\n" % (ii, jj, k * 0.5, r0))
    #LTfile.write("   bond_coeff  @bond:%d_%d  %f   %f\n" % (ii, jj, k * 0.5, r0))
np.savetxt("K_const_ind.txt",k_const_ind,fmt='%d')
LTfile.write(" }\n")

LTfile.write("}\n")
LTfile.write("pep = new CG[%d]"%(n_chains)) #TODO: define 'n_chains' above
