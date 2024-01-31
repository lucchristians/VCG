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

def create_hist(data,n_frames,Np):
    shape=np.shape(data)[0]
    rlo=0.0
    rhi=20.0
    g , x = np.histogram(data,bins=1000,range=(rlo,rhi))
    #print(np.sum(g))
    return x[:-1], g/n_frames/Np

def pair_hist(txt, pos_array,atomtype):
    ranges=[]
    probs=[] #TODO: try to speed up
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
                #bc=b[:,j,:]-c[:,k,:] #TODO: changed due to limitations with model set up
                d.append(b[:,j,:]-c[:,k,:])
                #d.append(bc-L*np.rint(bc/L))
                Np+=1
        d=np.swapaxes(np.array(d),0,1) #shape = [n_frames,Np,3]
        
        
        avg_mode=False
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
    ranges=np.array(ranges)
    probs=np.array(probs)
    np.savetxt("rcut_stacked/CG_hist.txt",probs) 
    np.savetxt("rcut_stacked/CG_ranges.txt",ranges) 

    return ranges, probs
#determing vs-connected pairs
vs_loc=[60]
#data=np.loadtxt("rcut_stacked/result.txt",skiprows=1) #TODO: representthis differently
data=[]
array_i=np.arange(1,22,dtype=int)
array_j=np.array([53])
for ind_i in array_i:
    for ind_j in array_j:
        data.append([ind_i,ind_j])
data=np.array(data)        
#inds=data[:,:2]
new_inds=data
#new_inds=[]

#r0=[] #TODO: populate
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


if(not os.path.exists("rcut_stacked/CG_ranges.txt") and not os.path.exists("rcut_stacked/CG_hist.txt")): 
    traj=Trajectory("/mnt/data0/decompose_forces/Qwt/run0/AA_stack_retyped.lammpstrj", fmt='lammpstrj') 
    pos_array=[]
    while(traj.read_frame()):
        atomtype=traj.t
        pos=traj.x
        pos_array.append(np.copy(pos))
    n_frames=len(pos_array)
    pos_array=np.array(pos_array)
    print(atomtype)
    print(np.shape(atomtype))
    ranges, probs = pair_hist(new_inds,pos_array,atomtype)
else:
    ranges=np.loadtxt("rcut_stacked/CG_ranges.txt") 
    probs =np.loadtxt("rcut_stacked/CG_hist.txt") 

#curve fit to find r0's and k's
constants=[]
#x=np.array([ranges,r0])

for i in range(len(new_inds)):
    
    inits=[0.005,10.0]
    print(new_inds[i])
    hist_noise_reduced= sig.savgol_filter(probs[i],51,5)   
    #plt.plot(probs[i],linewidth=1)
    #plt.plot(hist_noise_reduced,linewidth=1)
    #plt.show()
    #exit()
    peaks, properties = sig.find_peaks(hist_noise_reduced,prominence=0.000005) #TODO: determine if prominence should be considered on a cause by case basis
    #print(peaks)
    max_peak=-1
    max_peak_ind=-1
    num=0
    input_mode=False
    if(input_mode):
        for peak in peaks:
            print("peak %d | %d @ %7.4f A with height %7.4f"%(num,peak,ranges[i,peak],hist_noise_reduced[peak]))
        peakx=float(input("what peak position will be used: "))
        #TODO: ignore if negative
        if(peakx < 0):
            continue
    #TODO: get ind if not
        max_peak_ind=len(ranges[i,peakx>ranges[i,:]])
    else:    
        for peak in peaks:
            height=hist_noise_reduced[peak]
            if(max_peak < height): #TODO: consider how best to represent this
                max_peak=height
                max_peak_ind=peak
            #plt.plot(x[peak+1],hist_noise_reduced[peak],'r*') print(peaks)
    
    peak_ind=max_peak_ind
    peakx=ranges[i,peak_ind]
    #peak_ind=peaks[0] #TODO: this doesn't seem realistic
    #peakx=ranges[i,peak_ind]
    x=np.array([ranges[i],peakx*np.ones(len(ranges[i]))])
    #print(peak_ind)
    #print(peakx)
    #peakx=ranges[i,peak_ind]
    start=peak_ind-30
    end=peak_ind+30
    opt, _ =op.curve_fit(harmonic,x[:,start:end],probs[i,start:end],p0=inits)
    plt.figure()
    plt.plot(ranges[i],hist_noise_reduced)
    plt.plot(ranges[i],harmonic(x,opt[0],opt[1]))
    plt.title(i)
    plt.show()
    if(peakx<15):
        constants.append([new_inds[i,0],60,opt[0],peakx]) #TODO: filter peaks that make no sense; determine good criteria
constants=np.array(constants)
np.savetxt("rcut_stacked/constants.txt",constants) #TODO: produce AA_hist file for this virtual site
