import numpy  as np
import scipy.optimize as op
import os
import scipy.signal as sig
import matplotlib.pyplot as plt
from mscg import *
L=255.0

def harmonic(x,k):
    r=x[0]
    r0=x[1]
    return np.exp(-k*(r-r0)**2)

def create_hist(data,n_frames,Np):
    shape=np.shape(data)[0]
    rlo=0.0
    rhi=20.0
    g , x = np.histogram(data,bins=1000,range=(rlo,rhi))
    #print(np.sum(g))
    return x[:-1], g/n_frames/Np

def pair_hist(txt, pos_array):
    ranges=[]
    probs=[] #TODO: try to speed up
    for i in range(len(txt)):
        # Shape array from [frames, number of types, 3] to [frames * number of types, 3]
        atom1=txt[i,0]
        atom2=txt[i,1]
        b = pos_array[:,atomtype==atom1,:]
        c = pos_array[:,atomtype==atom2,:]
        d = []

        # Isolate interaction interface associated interactions
        Np=0
        for j in range(np.shape(b)[1]):
            for k in range(np.shape(c)[1]):
                if(not j==k):
                    continue #TODO: min img conv
                #if(j==(k+1)%np.shape(b)[1]): # Compares chain to adjacent vs-interacting chain
                bc=b[:,j,:]-c[:,k,:]
                d.append(bc-L*np.rint(bc/L))
                Np+=1
        d=np.swapaxes(np.array(d),0,1)
        e = np.linalg.norm(d,axis=2).flatten()
        f = e
        #print(np.shape(f))
        x, h = create_hist(f,n_frames,Np)
        ranges.append(x)
        probs.append(np.copy(h))
    ranges=np.array(ranges)
    probs=np.array(probs)
    np.savetxt("rcut_combine/CG_hist.txt",probs)
    np.savetxt("rcut_combine/CG_ranges.txt",ranges)

    return ranges, probs
#determing vs-connected pairs
vs_loc=[55,56,57,58,59]
data=np.loadtxt("rcut_combine/result.txt",skiprows=1)
inds=data[:,:2]
new_inds=[]
r0=[]
for i in range(len(inds)):
    vs_present=False
    for j in range(len(vs_loc)):
        if(inds[i,0]==vs_loc[j] or inds[i,1]==vs_loc[j]):
            vs_present=True
            break
    if(vs_present):
        #r0.append(data[i,2]*np.ones(1000))
        new_inds.append(inds[i])
new_inds=np.array(new_inds)

#generate hists if they dont exist/ load them in if they dont
traj=Trajectory("vs_m.lammpstrj", fmt='lammpstrj')
pos_array=[]
while(traj.read_frame()):
    atomtype=traj.t
    pos=traj.x
    pos_array.append(np.copy(pos))
n_frames=len(pos_array)
pos_array=np.array(pos_array)
print(n_frames)
if(not os.path.exists("rcut_combine/CG_ranges.txt") and not os.path.exists("rcut_combine/CG_hist.txt")):
    ranges, probs = pair_hist(new_inds,pos_array)
else:
    ranges=np.loadtxt("rcut_combine/CG_ranges.txt")
    probs =np.loadtxt("rcut_combine/CG_hist.txt")

#curve fit to find r0's and k's
constants=[]
#x=np.array([ranges,r0])

for i in range(len(new_inds)):
    inits=[0.005,10.0]
    #TODO: limit range to gaussian
    print(new_inds[i])
    hist_noise_reduced= sig.savgol_filter(probs[i],51,5)   
    #plt.plot(probs[i],linewidth=1)
    #plt.plot(hist_noise_reduced,linewidth=1)
    #plt.show()
    #exit()
    peaks, properties = sig.find_peaks(hist_noise_reduced,prominence=0.005)
    max_peak=0
    for peak in peaks:
        height=hist_noise_reduced[peak]
        if(max_peak < height):
            max_peak=height
            max_peak_ind=peak
        #plt.plot(x[peak+1],hist_noise_reduced[peak],'r*') print(peaks)
    
    peak_ind=max_peak_ind
    
    peakx=ranges[i,peak_ind]
    x=np.array([ranges[i],peakx*np.ones(len(ranges[i]))])
    #r0.append(data[i,2]*np.ones(1000))
    start=peak_ind-30
    end=peak_ind+30
    opt, _ =op.curve_fit(harmonic,x[:,start:end],probs[i,start:end],p0=inits)
    constants.append([new_inds[i,0],new_inds[i,1],opt[0],peakx,opt[1]])
constants=np.array(constants)
np.savetxt("rcut_combine/constants.txt",constants)
