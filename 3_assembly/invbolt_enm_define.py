import numpy  as np
import scipy.optimize as op
import os
import scipy.signal as sig
import matplotlib.pyplot as plt
from mscg import *
#L=255.0
kb= 0.0019872041 #[units]
T=277 #K
beta=1/(kb*T)

def harmonic(x,k,norm):
    r=x[0]
    r0=x[1]
    return norm*np.exp(-beta*k*(r-r0)**2)

def find_frames_in_range(pdist,min_,max_):
    temp=[]
    for i in range(len(pdist)):
        if(pdist[i]>min_ and pdist[i]<max_):
            temp.append(i)
    return np.array(temp)

def create_hist(data,n_frames,Np):
    shape=np.shape(data)[0]
    rlo=0.0
    rhi=20.0
    g , x = np.histogram(data,bins=1000,range=(rlo,rhi))
    #print(np.sum(g))
    return x[:-1], g/n_frames/Np

def pair_hist(txt, pos_array,atomtype):
    ranges=[]
    probs=[]
    f_list=[] 
    for i in range(len(txt)):
        # Shape array from [frames, number of types, 3] to [frames * number of types, 3]
        atom1=txt[i,0]
        atom2=txt[i,1]
        b = pos_array[:,atomtype==atom1,:]
        c = pos_array[:,atomtype==atom2,:]
        d = []
        #print(np.shape(b)[1])
        #print(np.shape(c)[1])
        # Isolate interaction interface associated interactions
        Np=0
        for j in range(np.shape(b)[1]):
            for k in range(np.shape(c)[1]):
                if(j==k): #are there ways to make this realistic
                    continue 
                #if(j==(k+1)%np.shape(b)[1]): # Compares chain to adjacent vs-interacting chain
                #bc=b[:,j,:]-c[:,k,:] 
                d.append(b[:,j,:]-c[:,k,:])
                #d.append(bc-L*np.rint(bc/L))
                Np+=1
        d=np.swapaxes(np.array(d),0,1) #shape = [n_frames,Np,3]
        
        
        avg_mode=True
        if(avg_mode):
            e = np.linalg.norm(d,axis=2).flatten()
            f=e
        else:
            e = np.linalg.norm(d,axis=2)
            f = np.min(e,axis=1)
        
        #Np=10
        #print(np.shape(f))
        x, h = create_hist(f,n_frames,Np)
        ranges.append(x)
        probs.append(np.copy(h))
        f_list.append(np.copy(f))
    ranges=np.array(ranges)
    probs=np.array(probs)
    np.savetxt("rcut_stacked/CG_pdist.txt",f_list)
    np.savetxt("rcut_stacked/CG_hist.txt",probs) 
    np.savetxt("rcut_stacked/CG_ranges.txt",ranges) 

    return ranges, probs, f_list #TODO: change instanences of this to include added out
#determing vs-connected pairs
output_ind={1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10,12:12,13:13,53:60,54:61}
#data=np.loadtxt("rcut_stacked/result.txt",skiprows=1) 
data=[]
#array_i=np.array([1,2,3,4,5,6,7,8,9,10])
array_i=np.arange(28,59)
array_j=np.array([62,63,64])
vs_pairs=np.array([[2,53],[1,54]])
for ind_i in array_i:
    for ind_j in array_j:
        if(ind_i==ind_j): #prevents overlap
            continue
        data.append([ind_i,ind_j])
data=np.array(data)        
#inds=data[:,:2]
new_inds=data
#new_inds=[]

#r0=[] 
'''
for i in range(len(inds)):
    vs_present=False
    for j in range(len(vs_loc)):
        if(inds[i,0]==vs_loc[j] or inds[i,1]==vs_loc[j]):
            vs_present=True
            break
    if(vs_present):
        r0.append(data[i,2]*np.ones(1000)) #TODO: parameterize this based on peak location
        new_inds.append(inds[i])
new_inds=np.array(new_inds)
'''
#generate hists if they dont exist/ load them in if they dont

'''
if(not os.path.exists("CG_ranges.txt") and not os.path.exists("CG_hist.txt")): 
    traj=Trajectory("AA.lammpstrj", fmt='lammpstrj')
    pos_array=[]
    while(traj.read_frame()):
        atomtype=traj.t
        pos=traj.x
        pos_array.append(np.copy(pos))
    n_frames=len(pos_array)
    pos_array=np.array(pos_array)[20000:] #ignore equilibration time
    print(atomtype)
    print(np.shape(atomtype))
    np.savetxt("CG_atomt.txt",atomtype)
    ranges, probs, xyz = pair_hist(new_inds,pos_array,atomtype)
else:
    ranges=np.loadtxt("CG_ranges.txt") 
    probs =np.loadtxt("CG_hist.txt") 
    xyz =np.loadtxt("CG_pdist.txt")
    atomtype =np.loadtxt("CG_atomt.txt")
'''    
traj=Trajectory("hp_vs_repl.lammpstrj", fmt='lammpstrj')
pos_array=[]
while(traj.read_frame()):
    atomtype=traj.t
    pos=traj.x
    pos_array.append(np.copy(pos))
n_frames=len(pos_array)
pos_array=np.array(pos_array)#[20000:] #ignore equilibration time
#curve fit to find r0's and k's
constants=[]
#x=np.array([ranges,r0])
redo_inds=[]
peak_finding_data=[]
for i in range(len(new_inds)):
    # Shape array from [frames, number of types, 3] to [frames * number of types, 3]
    atom1=new_inds[i,0]
    atom2=new_inds[i,1]
    b = pos_array[:,atomtype==atom1,:]
    c = pos_array[:,atomtype==atom2,:]
    #print(b.shape)
    #print(c.shape)
    d = []
    # Isolate interaction interface associated interactions
    Np=0
    for j in range(np.shape(b)[1]):
        for k in range(np.shape(c)[1]):
            if(j==k): #are there ways to make this realistic
                #continue 
                d.append(b[:,j,:]-c[:,k,:])
    d=np.swapaxes(np.array(d),0,1) #shape = [n_frames,Np,3]
    e = np.linalg.norm(d,axis=2)
    print(e)
    peakx=np.mean(e,axis=1)[0]
    print(peakx)
    #exit()
    if(peakx<7.4):
        constants.append([new_inds[i,0],new_inds[i,1],1.5,peakx]) #elastic network 


constants=np.array(constants)
np.savetxt("constants.txt",constants) 
