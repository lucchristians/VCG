import numpy as np
import mdtraj as md

# Input options
stagger=64.0
r=28.0
theta=50.0 /180*np.pi
alpha=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
anti_parallel = True
pdb2 = open("Qwt_v2.pdb",'w')
pdb_base = md.load("bot.pdb")
top_info = pdb_base.top.to_dataframe()
height=2
width=1
length=1
#height=5
#width=4
#length=4


# a rotation fucntion
def rotate(points, pitch, roll, yaw): #np.shape(points) = [n_atoms,3], others are constants
    cosa = np.cos(yaw)
    sina = np.sin(yaw)

    cosb = np.cos(pitch)
    sinb = np.sin(pitch)

    cosc = np.cos(roll)
    sinc = np.sin(roll)

    Axx = cosa*cosb
    Axy = cosa*sinb*sinc - sina*cosc
    Axz = cosa*sinb*cosc + sina*sinc

    Ayx = sina*cosb
    Ayy = sina*sinb*sinc + cosa*cosc
    Ayz = sina*sinb*cosc - cosa*sinc

    Azx = -sinb
    Azy = cosb*sinc
    Azz = cosb*cosc
    new_points = np.copy(points)
    for i in range(len(points)):
        px = points[i,0]
        py = points[i,1]
        pz = points[i,2]

        new_points[i,0] = Axx*px + Axy*py + Axz*pz
        new_points[i,1] = Ayx*px + Ayy*py + Ayz*pz
        new_points[i,2] = Azx*px + Azy*py + Azz*pz
    
    return new_points



# Interpret the coors of an input pdb
n_atoms = pdb_base.n_atoms
print(n_atoms)
n_chains = pdb_base.topology.n_chains
atoms_per_chain=int(n_atoms/n_chains)
print(atoms_per_chain)
xyz = pdb_base.xyz*10 # the original pdb was generated in terms of angstroms
xyz = xyz.squeeze()

# Center the pdb around the origin
shift = md.compute_center_of_mass(pdb_base)
print(shift)
xyz_centered = xyz -shift*10


#build the structure 2
index=1 # atom index
res=1
chain=0 #represents chain
for l in range(height): #height
    for k in range(width): #y dir
        for m in range(length): #z dir
            #TODO: find ways to fibrilizes for nonsquare lattice
            res=0
            
            comx=64.*l+500+(stagger*((k+m)%2))
            #print(np.pi*((k+m)%2))
            #print(comx)
            #TODO: determine effective radius of Qwt
            comz=r*m+500
            comy=-r*k+500
            if(anti_parallel):
                xyz_oriented =  rotate(xyz_centered,0.0,theta*((k+m)%2),np.pi*((k+m)%2)) 
            else:
                xyz_oriented =  rotate(xyz_centered,0.0,theta*((k+m)%2),0.0)
            prev_res_name=top_info[0]['resName'][0]
            for i in range(n_atoms):
                #res=1
                if(top_info[0]['name'][i]=='N'):
                    res+=1
                x=xyz_oriented[i,0] + comx
                y=xyz_oriented[i,1] + comy
                z=xyz_oriented[i,2] + comz
                pdb2.write("%-4s  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%("ATOM",index%100000,top_info[0]['name'][i],top_info[0]['resName'][i],alpha[chain%len(alpha)],res,x,y,z,0.00,0.00))
                index+=1
                prev_res_name=top_info[0]['name'][i]
                
                
                if(i%atoms_per_chain == atoms_per_chain-1):
                    res=0
                    chain+=1
pdb2.write("END\n")
pdb2.close()

