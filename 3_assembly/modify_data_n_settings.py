import numpy as np
from mscg import *
import scipy.optimize as op
import sys
input_data=sys.argv[1]
output_data=sys.argv[2]
#find a position from 4 points and lengths
def fit(inpu,x,y,z): #shape of x: [4,3]
    pos=np.array([[x,y,z]])
    return np.linalg.norm(inpu-pos,axis=1)
    

def triangluate(xyz,lines):
    opt,_=op.curve_fit(fit,xyz,lines,p0=np.mean(xyz,axis=0))
    return opt

# read in location data


#read old data file
init=True
masses_section=False
atoms_section=False
bonds_section=False
charge=[]
masses=[]
pos=[]
bond_list=[]
with open(input_data,'r') as f:
    for line in f:
        
        dat=line.split()
        if(len(dat)==0):
            #print("")
            continue
        if(init and len(dat)>1):
            if(dat[1]=='atoms'):
                atoms=int(dat[0])
            elif(dat[1]=='bonds'):
                bonds=int(dat[0])
            elif(dat[1]=='angles'):
                angles=int(dat[0])
            elif(dat[1]=='dihedrals'):
                dihedrals=int(dat[0])
            elif(dat[1]=='impropers'):
                impropers=int(dat[0])
            elif(dat[1]=='atom'):
                atypes=int(dat[0])
            elif(dat[1]=='bond'):
                btypes=int(dat[0])
            elif(len(dat)>=3):
                if(dat[2]=='xlo'):
                    xhi=float(dat[1])
                if(dat[2]=='ylo'):
                    yhi=float(dat[1])
                if(dat[2]=='zlo'):
                    zhi=float(dat[1])

        if(dat[0]=='Masses'):
            init=False
            masses_section=True
        if(dat[0]=='Atoms'):
            masses_section=False
            atoms_section=True
            continue
        
        if(dat[0]=='Bonds'):
            atoms_section=False
            bonds_section=True
            continue
        if(masses_section and len(dat)>1):
            masses.append(float(dat[1]))
        elif(atoms_section and len(dat)>6):
            charge.append([dat[0],dat[1],dat[2],dat[3]])
            pos.append([float(dat[4]),float(dat[5]),float(dat[6])])
        elif(bonds_section and len(dat)>3):
            bond_list.append([int(dat[0]),int(dat[1]),int(dat[2]),int(dat[3])])
        #else:
            #print(line[:-1])
masses=np.array(masses).astype(float)
charge=np.array(charge).astype(float)
pos=np.array(pos)
bond_list=np.array(bond_list)

#get constants
consts=np.loadtxt('constants.txt')
#create separate settings file based on information
init_bond=btypes+1
new_bonds=[]
ff=open('hydrophobic.ff','w')
#determine new virtual sites
unique=[]
for con in consts:
    if int(con[1]) not in unique:
        unique.append(int(con[1]))

print(unique)
#exit()
soft_repl=np.loadtxt('soft_const.txt')

for i, u in enumerate(unique): #TODO: 
    ff.write('pair_coeff * %d soft 0 5.0\n'%(soft_repl[i,0]))
    ff.write('pair_coeff %d %d soft %6.3f %6.3f\n'%(int(soft_repl[i,0]),int(soft_repl[i,1]),soft_repl[i,2],soft_repl[i,3]))
ff.write("\n")
for i in range(len(consts)):
    non=init_bond+i
    new_bonds.append([non,consts[i,0],consts[i,1]])
    ff.write('bond_coeff %d %7.4f %7.4f\n'%(non,consts[i,2],consts[i,3]))
new_bonds=np.array(new_bonds)
ff.close()

#triangulate positions of added vsites
locations=[]
for u in unique:
    data=consts[u==consts[:,1].astype(int),:]
    points=[]
    for i in range(len(data)):
        points.append(pos[charge[:,2].astype(int)==data[i,0].astype(int),:])
    points=np.array(points)
    print(points.shape)
    lines=data[:,3]
    loca=[]
    for i in range(points.shape[1]):
        loca.append(triangluate(points[:,i,:],lines))
    locations.append(loca)
locations=np.array(locations)
print(locations.shape)
#exit()

#add locations to pos,bonds, charge arrays
chains=int(len(charge)/atypes)
new_charge=[]
new_pos=[]
new_bond_list=[]
mono_charge=charge[:atypes]
#mono_pos=pos[:atoms]
mono_bonds=bond_list[:btypes]
for i in range(chains):
    print(i)
    #copy monomer
    temp_charge=np.copy(mono_charge)
    temp_pos=np.copy(pos[i*atypes:(i+1)*atypes])
    temp_bond=np.copy(mono_bonds)
    #manipulate
    temp_charge[:,0]=temp_charge[:,0]+i*(atypes+len(locations))
    temp_charge[:,1]=temp_charge[:,1]+i
    temp_bond[:,0]=temp_bond[:,0]+i*(btypes+len(consts))
    temp_bond[:,2]=temp_bond[:,2]+i*(atypes+len(locations))
    temp_bond[:,3]=temp_bond[:,3]+i*(atypes+len(locations))

    charge_loc_temp=[]
    pos_loc_temp=np.copy(locations[:,i,:])
    for j in range(len(locations)):
        charge_loc_temp.append([float(unique[j])+i*(atypes+len(locations)),1+i,unique[j],0.0000])
    bond_loc_temp=[]
    for j in range(len(consts)):
        bond_loc_temp.append([new_bonds[j,0]+btypes*i,new_bonds[j,0],new_bonds[j,1]+(atypes+len(locations))*i,new_bonds[j,2]+(atypes+len(locations))*i])

    charge_loc_temp=np.array(charge_loc_temp)
    temp_charge=np.append(temp_charge,charge_loc_temp,axis=0)
    temp_pos=np.append(temp_pos,pos_loc_temp,axis=0)
    temp_bond=np.append(temp_bond,bond_loc_temp,axis=0)
    
    new_charge.append(np.copy(temp_charge))
    new_pos.append(np.copy(temp_pos))
    new_bond_list.append(np.copy(temp_bond))
new_charge=np.vstack(new_charge)
new_pos=np.vstack(new_pos)
new_bond_list=np.vstack(new_bond_list)
atoms=new_charge.shape[0]
#new_pos.shape)
bonds=new_bond_list.shape[0]

#exit()

#create new data file with added vsites

ndf=open(output_data,'w')

ndf.write('LAMMPS Description\n')
ndf.write('\n')
ndf.write('%d atoms\n'%(int(atoms)))
ndf.write('%d bonds\n'%(int(bonds)))
ndf.write('%d angles\n'%(int(angles)))
ndf.write('%d dihedrals\n'%(int(dihedrals)))
ndf.write('%d impropers\n'%(int(impropers)))
ndf.write('\n')
ndf.write('%d atom types\n'%(int(atypes+len(unique))))
ndf.write('%d bond types\n'%(int(btypes+len(consts))))
ndf.write('\n')
ndf.write('%7.4f %7.4f xlo xhi\n'%(0,xhi))
ndf.write('%7.4f %7.4f ylo yhi\n'%(0,yhi))
ndf.write('%7.4f %7.4f zlo zhi\n'%(0,zhi))
ndf.write('\n')
ndf.write('Masses\n')
ndf.write('\n')
for i in range(len(masses)):
    ndf.write('%d %10.6f # CG%d\n'%(i+1,float(masses[i]),i+1))
for i in range(len(unique)):
    ndf.write('%d %10.6f # CG%d\n'%(int(unique[i]),50,int(unique[i])))
ndf.write('\n')
ndf.write('Atoms\n')
ndf.write('\n')
for i in range(len(new_pos)):
    ndf.write('%d %d %d %8.6f %6.3f %6.3f %6.3f\n'%(int(new_charge[i,0]),int(new_charge[i,1]),int(new_charge[i,2]),new_charge[i,3],new_pos[i,0],new_pos[i,1],new_pos[i,2]))
ndf.write('\n')
ndf.write('Bonds\n')
ndf.write('\n')
for i in range(len(new_bond_list)):
    ndf.write('%d %d %d %d\n'%(new_bond_list[i,0],new_bond_list[i,1],new_bond_list[i,2],new_bond_list[i,3]))
ndf.close()
