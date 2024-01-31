import numpy as np
import matplotlib.pyplot as plt
import sys
import os
print("start")
binwidth=0.05 #A
max=20 #20.0
bins=int(max/binwidth)
colors=['#CCCCCC','#E83274','#4FA3EC']

def compile_data(path,pid,eid):
    data=[]
    for k in range(5):
        data_t=[]
        for i in range(4):
            if(i==1):
                continue
            for j in range(pid,eid):
            
                #print(i,j,k)
                data_t.append(np.load(f'{path}{i+1}/pair{k+1}_step{j}.npy'))
        data.append(np.hstack(data_t))
    print(len(data))
    return data

def create_hist(data,lens,Norm):
    #shape=np.shape(data)[0]
    rlo=0.0
    rhi=max
    xs=[]
    hists=[]
    for i in range(lens):
        g , x = np.histogram(data[i],bins=bins,range=(rlo,rhi))
        xs.append(np.copy((x[:-1]+x[1:])/2.0))
        hists.append(np.copy(g/np.sum(g)*norms[i]))
    return xs, hists 

probs_AA = np.loadtxt("AA_hist.txt")
print(np.shape(probs_AA))
if(len(np.shape(probs_AA))==1):
    bins=len(probs_AA)
    probs_AA = probs_AA.reshape(1,bins)
else:
    bins=len(probs_AA[0])

ratio=1/5
inch=16
size=[inch,inch*ratio]
#fig = plt.figure(figsize=size)
print("formatting")
plt.style.use("~/fig_format.mplstyle") #basic formatting
fig, axs = plt.subplots(nrows=1,ncols=5,sharey=True,figsize=size)
plt.subplots_adjust(wspace=0.01)

x_=np.arange(0.05,20.0001,0.05)
norms=[]
for j in range(len(probs_AA)-2):
                #plt.figure(j)
                axs[j].plot(x_,probs_AA[j],'-',color=colors[0],label="AA") #fill
                norms.append(np.sum(probs_AA[j]))
                #ax = plt.axes()
                d=np.zeros(len(probs_AA[j]))
                axs[j].fill_between(x_,probs_AA[j], where=probs_AA[j]>=d, interpolate=True, color=colors[0])
                

print("loading CG")
vs_rep=['novs','vs']
for vs in vs_rep:
    
            if(vs=='vs'): #vsite
                path='vs_intra_repl_'
                i=2 #color
                txt='VCG'
            else: #gauss
                path='novs_intra_pak_'
                i=1 #color
                txt='CG'
            data=compile_data(path,4,11)
            data_len=len(data)
            x, temp = create_hist(data,data_len,norms)    
            #x, temp = set_average([num],path)
            for j in range(len(probs_AA)-2):
                #plt.figure(j)
                axs[j].plot(x[j],np.copy(temp[j]),'-',color=colors[i],linewidth=1,label=txt)
                axs[j].text(3.41,0.011,'VS%d'%(j+1),fontsize=12)
print("does it get here")
for j in range(len(probs_AA)-2):
    #plt.figure(j)
    if(j==2):
        axs[j].set_xlabel("Distance ($\mathrm{\AA}$)",fontsize=18,fontname='Arial')
    if(j==0):
        axs[j].set_ylabel("Probability",fontsize=18,labelpad=10)
    if(j==4):
        axs[j].legend()
    #plt.figure(j).subplots_adjust(left=0.16)
    axs[j].set_xlim(3,17)
    axs[j].set_ylim(0.0,0.0125)
#plt.tight_layout()
plt.subplots_adjust(bottom=0.27)
plt.savefig("Assemb_hist.png")
