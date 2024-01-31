import numpy as np
import sys
groPath=sys.argv[1]

# load gro file
with open(groPath) as f:
    count=1
    atoms=0
    xyz=[]
    inds=[]
    for line in f:
        if(count>2):
            if(line[0:3]=='END'):
                break
            elif(line[0:3]=='TER'):
                continue
            else:
                atomtype=line[13:15]
                if(atomtype=='CA'):
                    #print(atomtype)
                    ind=int(line[15:20])
                    x=float(line[20:28])
                    y=float(line[28:36])
                    z=float(line[36:44])
                    inds.append(ind)
                    xyz.append([x,y,z])
        count+=1
    xyz=np.array(xyz)
    inds=np.array(inds)
    atoms=len(xyz)
#determine length of a coiled-coil
'''
length_CC=int(atoms/2) #TODO: make work for 18-36 and 36-54

#print(np.shape(xyz))
#determine center of masses for both chains
CC1_xyz=xyz[:length_CC,:]
CC1_ind=inds[:length_CC]
CC2_xyz=xyz[length_CC:,:]
CC2_ind=inds[length_CC:]
'''
length_CC_top=[]
length_CC_bottom=[]
for i in range(5):
    length_CC_top.append(np.arange(18+54*i,36+54*i))
    length_CC_bottom.append(np.arange(36+54*i,54+54*i))
length_CC_top=np.array(length_CC_top)
length_CC_bottom=np.array(length_CC_bottom)    
#len_CC_top=int(len(length_CC_top)/2)
#len_CC_bottom=int(len(length_CC_bottom)/2)
len_CC_top=[0,1,3,4]
len_CC_bottom=[2]
#CC1_xyz_top=xyz[length_CC_top[:len_CC_top],:]
#CC1_ind_top=inds[length_CC_top[:len_CC_top]]
#CC2_xyz_top=xyz[length_CC_top[len_CC_top:],:]
#CC2_ind_top=inds[length_CC_top[len_CC_top:]]
CC1_xyz_top=xyz[length_CC_top[len_CC_bottom].flatten(),:]
CC1_ind_top=inds[length_CC_top[len_CC_bottom].flatten()]
CC2_xyz_top=xyz[length_CC_top[len_CC_top].flatten(),:]
CC2_ind_top=inds[length_CC_top[len_CC_top].flatten()]
#print(CC1_ind_top)
#print(CC2_ind_top)
#CC1_xyz_bottom=xyz[length_CC_bottom[:len_CC_bottom],:]
#CC1_ind_bottom=inds[length_CC_bottom[:len_CC_bottom]]
#CC2_xyz_bottom=xyz[length_CC_bottom[len_CC_bottom:],:]
#CC2_ind_bottom=inds[length_CC_bottom[len_CC_bottom:]]
CC1_xyz_bottom=xyz[length_CC_bottom[len_CC_bottom].flatten(),:]
CC1_ind_bottom=inds[length_CC_bottom[len_CC_bottom].flatten()]
CC2_xyz_bottom=xyz[length_CC_bottom[len_CC_top].flatten(),:]
CC2_ind_bottom=inds[length_CC_bottom[len_CC_top].flatten()]


CC1_com_top=np.sum(CC1_xyz_top,axis=0)/len(CC1_xyz_top)
CC2_com_top=np.sum(CC2_xyz_top,axis=0)/len(CC2_xyz_top)
direction_top=CC1_com_top-CC2_com_top

CC1_com_bottom=np.sum(CC1_xyz_bottom,axis=0)/len(CC1_xyz_bottom)
CC2_com_bottom=np.sum(CC2_xyz_bottom,axis=0)/len(CC2_xyz_bottom)
direction_bottom=CC1_com_bottom-CC2_com_bottom

#print(np.shape(CC1_xyz))
#print(np.shape(CC2_xyz))

#determine atom distance 
rij_min_top=9999
i_min_top=-1
j_min_top=-1
for i in range(len(CC1_xyz_top)):
    for j in range(len(CC2_xyz_top)):
        dxyz_top=CC1_xyz_top[i]-CC2_xyz_top[j]
        if(np.any((dxyz_top/direction_top)<0.0)):
            continue
        rij_top=np.linalg.norm(dxyz_top)
        if(rij_top<rij_min_top):
            #print(rij)
            rij_min_top=rij_top
            i_min_top=CC1_ind_top[i]
            j_min_top=CC2_ind_top[j]

rij_min_bottom=9999
i_min_bottom=-1
j_min_bottom=-1
for i in range(len(CC1_xyz_bottom)):
    for j in range(len(CC2_xyz_bottom)):
        dxyz_bottom=CC1_xyz_bottom[i]-CC2_xyz_bottom[j]
        if(np.any((dxyz_bottom/direction_bottom)<0.0)):
            continue
        rij_bottom=np.linalg.norm(dxyz_bottom)
        if(rij_bottom<rij_min_bottom):
            #print(rij)
            rij_min_bottom=rij_bottom
            i_min_bottom=CC1_ind_bottom[i]
            j_min_bottom=CC2_ind_bottom[j]


print(i_min_top,j_min_top,i_min_bottom,j_min_bottom) #TODO: what script do these numbers go with        

