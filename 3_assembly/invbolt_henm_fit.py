import numpy  as np
import scipy.optimize as op
import os
import scipy.signal as sig
import matplotlib.pyplot as plt
from mscg import *
#L=255.0
kb= 0.0019872041 #[kcal/mol/K]
T=277 #K
beta=1/(kb*T)

def harmonic(x,k,norm):
    r=x[0]
    r0=x[1]
    return norm*np.exp(-beta*k*(r-r0)**2)

def soft(x,D,norm):
    r=x[0]
    rc=x[1]
    return norm*np.exp(-beta*D*(1+np.cos(np.pi*r/rc)))
'''
def soft(x,D,rc,norm):
    return norm*np.exp(-beta*D*(1+np.cos(np.pi*x/rc)))
'''
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

#determing vs-connected pairs
#output_ind={1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10,12:12,13:13,53:60,54:61}
mapped={55:62,56:63,57:64,58:65}
#data=np.loadtxt("rcut_stacked/result.txt",skiprows=1) 
data=[]
#array_i=np.array([1,2,3,4,5,6,7,8,9,10])
array_i=np.arange(28,54)
array_j=np.array([55,56,57,58])
vs_pairs=np.array([[2,53],[1,54]])
for ind_i in array_i:
    for ind_j in array_j:
        if(ind_i==ind_j): #prevents overlap
            continue
        data.append([ind_i,ind_j])
data=np.array(data)        
#inds=data[:,:2]
new_inds=data

#load in trajectory
traj=Trajectory("AA_repl.lammpstrj", fmt='lammpstrj') #repl-vs in CG-mapped atomistic dataset
pos_array=[]
while(traj.read_frame()):
    atomtype=traj.t
    pos=traj.x
    pos_array.append(np.copy(pos))
n_frames=len(pos_array)
pos_array=np.array(pos_array)#[20000:] #ignore equilibration time

#curve fit to find r0's and k's for bonded potentials
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
                Np+=1
    d=np.swapaxes(np.array(d),0,1) #shape = [n_frames,Np,3]
    e = np.linalg.norm(d,axis=2).flatten()
    f=e
    #print(e)
    x_new, hist_new = create_hist(f,n_frames,Np)
    #else:
    #    x_new, hist_new = create_hist(pdists[ind],n_frames,Np)
    
    hist_noise_reduced= sig.savgol_filter(hist_new,51,5)
        
    peaks, properties = sig.find_peaks(hist_noise_reduced,prominence=0.000005)
    max_peak=-1
    max_peak_ind=-1
    for peak in peaks:
        height=hist_noise_reduced[peak]
        if(max_peak < height): #TODO: consider how best to represent this
            max_peak=height
            max_peak_ind=peak
        
    if(max_peak_ind<30):
        peak_ind=30
    else:
        peak_ind=max_peak_ind
    #peakx=np.mean(e,axis=1)[0]
    peakx=x_new[peak_ind]
    x=np.array([x_new,peakx*np.ones(len(x_new))])
    start=peak_ind-30
    end=peak_ind+30
        
    #fit model based based on peak   
    inits=[0.005,10.0]
    opt, _ =op.curve_fit(harmonic,x[:,start:end],hist_new[start:end],p0=inits)
    print(peakx)
    #exit()
    if(new_inds[i,1]>54):
        out_ind=mapped[new_inds[i,1]]
    else:
        out_ind=new_inds[i,1]
    
    if(peakx<10.0 and opt[0]>0.1):
        plt.figure()
        plt.plot(x_new,hist_new)
        plt.plot(x_new,harmonic(x,opt[0],opt[1]))
        plt.title("%d_%d"%(new_inds[i,0],new_inds[i,1]))
        plt.savefig("pair%d_%d.png"%(new_inds[i,0],new_inds[i,1]))
        plt.close()
        constants.append([new_inds[i,0],np.copy(out_ind),opt[0],peakx]) #elastic network 
np.savetxt("constants.txt",constants,fmt='%6.3f') 
#determine soft exclusions for 5-mer
#curve fit to find r0's and k's for bonded potentials
soft_const=[]
#x=np.array([ranges,r0])
peak_finding_data=[]
repl_pairs_ind=np.array([[55,55],[56,56],[57,57],[58,58]])
for i in range(len(repl_pairs_ind)):
    # Shape array from [frames, number of types, 3] to [frames * number of types, 3]
    atom1=repl_pairs_ind[i,0]
    atom2=repl_pairs_ind[i,1]
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
                continue 
            d.append(b[:,j,:]-c[:,k,:])
            Np+=1
    d=np.swapaxes(np.array(d),0,1) #shape = [n_frames,Np,3]
    e = np.linalg.norm(d,axis=2).flatten()
    f=e
    #print(e)
    x_new, hist_new = create_hist(f,n_frames,Np)
    #else:
    #    x_new, hist_new = create_hist(pdists[ind],n_frames,Np)
    
    hist_noise_reduced= sig.savgol_filter(hist_new,51,5)
        
    peaks, properties = sig.find_peaks(hist_noise_reduced,prominence=0.000005)
    max_peak=-1
    max_peak_ind=-1
    for peak in peaks:
        height=hist_noise_reduced[peak]
        if(max_peak < height): #TODO: consider how best to represent this
            max_peak=height
            max_peak_ind=peak
        
    if(max_peak_ind<30):
        peak_ind=30
    else:
        peak_ind=max_peak_ind
    #peakx=np.mean(e,axis=1)[0]
    peakx=x_new[peak_ind]
    x=np.array([x_new,peakx*np.ones(len(x_new))])
    #start=peak_ind-30
    #end=peak_ind+30
    start=0
    end=peak_ind    
    #fit model based based on peak   
    #inits=[0.008,peakx,5.0]
    #opt, _ =op.curve_fit(soft,x_new[start:end],hist_new[start:end],p0=inits)
    inits=[0.005,10.0]
    x=np.array([x_new,peakx*np.ones(len(x_new))])
    opt, _ =op.curve_fit(soft,x[:,start:end],hist_new[start:end],p0=inits)
    print(peakx)
    #exit()
    out_ind=np.zeros(2)
    for j in range(2):
        if(repl_pairs_ind[i,j]>54):
            out_ind[j]=mapped[repl_pairs_ind[i,j]]
        else:
            out_ind[j]=repl_pairs_ind[i,j]
    
    #if(peakx<10.0 and opt[0]>0.1):
    plt.figure()
    plt.plot(x_new,hist_new)
    #plt.plot(x_new,soft(x_new,opt[0],opt[1],opt[2]))
    plt.plot(x_new,soft(x,opt[0],opt[1]))
    plt.xlim(0,5)
    plt.title("%d_%d"%(repl_pairs_ind[i,0],repl_pairs_ind[i,1]))
    plt.savefig("soft%d_%d.png"%(repl_pairs_ind[i,0],repl_pairs_ind[i,1]))
    plt.close()
    soft_const.append([np.copy(out_ind[0]),np.copy(out_ind[1]),opt[0],peakx]) #elastic network 
soft_const=np.array(soft_const)
np.savetxt("soft_const.txt",soft_const,fmt='%6.3f') 
