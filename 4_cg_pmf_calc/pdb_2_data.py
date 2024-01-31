import numpy as np
import sys

pdb=sys.argv[1]
positions=[]

#getting charges
atoms_section=False
charge=[]
with open('outfile.data','r') as f:
    for line in f:
        
        dat=line.split()
        if(len(dat)==0):
            print("")
            continue
        if(dat[0]=='Atoms'):
            atoms_section=True
            continue
        
        if(dat[0]=='Bonds'):
            break
	
        if(atoms_section):
            charge.append([dat[0],dat[1],dat[2],dat[3]])
        else:
            print(line[:-1])
num=0
with open(pdb,'r') as f:
    for line in f:
        if(len(line)==0):
            continue

        if(line[:4]=='ATOM'):
            print(num+1,int(num/len(charge))+1,charge[num%len(charge)][2],charge[num%len(charge)][3],line[30:38],line[38:46],line[46:54])
            num+=1




