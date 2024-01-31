import numpy as np
import mdtraj as md
import sys
import matplotlib.pyplot as plt

target_traj=sys.argv[1]
target_top=sys.argv[2]
out_name=sys.argv[3]
def min_image(xyz,ref_ind,box):
    xyz_0=xyz[:,[ref_ind],:]
    #shift=np.array([[[60,60,60]]])-xyz0
    dr=xyz-xyz_0
    new_box=box.reshape([1,1,-1])
    dr = dr - new_box * np.rint( dr / new_box )
    return dr

def compute_theta(vector_1,vector_2):
    #print(vector_1,vector_2)
    #print(np.sum(vector_1*vector_2))
    #print(np.arccos((np.sum(vector_1*vector_2))/(np.linalg.norm(vector_1)*np.linalg.norm(vector_2))))
    #exit()
    return np.arccos((np.sum(vector_1*vector_2))/(np.linalg.norm(vector_1)*np.linalg.norm(vector_2)))

def sym_calc(xyz): #5-point input only
    #find center of geo
    com=np.mean(xyz,axis=0)
    #generate vectors from center of mass
    vectors=[]
    for points in xyz:
        vectors.append(points-com)
    vectors=np.array(vectors)
    #print(vectors)
    #print(vectors.shape)
    #exit()
    # find thetas between i and i+1%5
    sum_=np.array([0,0])
    for i in range(5):
        theta=compute_theta(vectors[i%5,:],vectors[(i+1)%5,:])
        #print(theta)
        sum_=sum_+np.array([np.cos(5*theta),np.sin(5*theta)])
        
    # sum it
    return sum_/5
    
def dist_calc(xyz): #5-point input only
    #find center of geo
    com=np.mean(xyz,axis=0)
    #generate vectors from center of mass
    vectors=[]
    for points in xyz:
        vectors.append(points-com)
    vectors=np.array(vectors)
    dist=np.linalg.norm(vectors,axis=-1)
    #print(np.shape(dist))
    return dist
    
    
targ=md.load(target_traj,top=target_top)
#refe=md.load(target_top)
ncg=int(targ.n_atoms/5)

sites=np.arange(18,49)
#sites=[37]
score=[]
xyz = min_image(targ.xyz,int(33+0*ncg),targ.unitcell_lengths[1]/2)
#xyz=targ.xyz
all_scores=[]
sym_s=[]
for site in sites:
    score=[]
    for f in range(targ.n_frames):
        s_inds=[]
        for i in range(0,5):
            s_inds.append(int(site+i*ncg))
        #print(s_inds)
        traj=targ.atom_slice(s_inds)    
        #score.append(md.acylindricity(traj)) 
        score.append(sym_calc(xyz[f,s_inds,:])[0]) 
        #temp=dist_calc(xyz[f,s_inds,:])
        #for t in temp:
        #    score.append(t) 
    all_scores.append(score)
    score=np.array(score)
    mean=np.mean(score,axis=0)
    std=np.std(score,axis=0)
    sym_s.append([mean,std])
sym_s=np.array(sym_s)
np.savetxt(out_name,[sym_s[:,0],sym_s[:,1]])
#plt.figure()
#plt.plot(sym_s[:,0,0])
#plt.savefig(out_name)
#for f in range(targ.n_frames):
'''
for site in sites:
    s_inds=[]
    for i in range(5,10):
        s_inds.append(int(site+i*ncg))
    traj=targ.atom_slice(s_inds)    
    score.append(md.acylindricity(traj)) 
    #score.append(sym_calc(xyz[f,s_inds,:])) 
#score=np.array(score)
score=np.hstack(score).flatten()
'''
all_=np.hstack(all_scores)
#print(np.shape(all_))
mean=np.mean(all_,axis=0)
std=np.std(all_,axis=0)
print("score: ",mean,"+/-",std) #,'+',mean[1],"+/-",std[1],"i")
#print("score: ",mean,"+/-",std,'nm') #,'+',mean[1],"+/-",std[1],"i")

#selection=targ.top.select('resid 18 to 48 and chainid 5 to 9')



#rmsd=md.rmsd(targ,refe,atom_indices=selection)

#print("RMSD: ",rmsd,"nm")

