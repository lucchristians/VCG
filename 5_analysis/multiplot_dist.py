import numpy as np
import matplotlib.pyplot as plt
import sys
import os

max=20.0
L=255.0
colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']
path_init = os.getcwd() 
consts=np.genfromtxt('reference_AA/const_AA.txt',ndmin=2)
probs_AA = np.loadtxt("reference_AA/AA_hist.txt")
#print(np.shape(probs_AA))
if(len(np.shape(probs_AA))==1):
    bins=len(probs_AA)
    probs_AA = probs_AA.reshape(1,bins)
else:
    bins=len(probs_AA[0])

def set_average(int_array, path):
    setn=[]
    #x=np.loadtxt(path+"/rem_iteration"+str(int_array[0])+"/CG_ranges.txt")
    x=np.loadtxt(path+"/rem_iteration"+str(int_array[0])+"/CG_r.txt")
    
    for i in int_array:
        #os.chdir(path+"/rem_iteration"+str(i))
        #traj=Trajectory(path+"/rem_iteration"+str(i)+"/sub1.lammpstrj", fmt='lammpstrj') #TODO: name is placeholder
        #x, temp=AA_2_CG(traj,probs_AA,txt,consts)
        temp=np.loadtxt(path+"/rem_iteration"+str(i)+"/CG_hist.txt")
        setn.append(np.copy(temp))
    setn=np.array(setn)

    setn_m=np.mean(setn,axis=0)
    return x, setn_m

vs_rep=['novs','vs']
intras=['','_bc']
dataset=['','_nointra']
#set_n=[]
x_=np.arange(0.05,20.0001,0.05)
for j in range(len(probs_AA)):
                plt.figure(j)
                plt.style.use("~/fig_format.mplstyle") #basic formatting
                plt.plot(x_,probs_AA[j],'-',color='#CCCCCC',label="AA") #fill
                ax = plt.axes()
                d=np.zeros(len(probs_AA[j]))
                ax.fill_between(x_,probs_AA[j], where=probs_AA[j]>=d, interpolate=True, color='#CCCCCC')
for vs in vs_rep:
    #for dat in dataset:
    #    for intra in intras:
            
            #for i in range(4):
            #path=path_init+"/rem_reps_"+vs+"_stacked"+dat+intra+"_"+str(i+1)
            if(vs=='vs'): #vsite
                path=path_init+"/vsite"
                i=1 #color
                num=243
                txt='VCG'
            else: #gauss
                path=path_init+'/gauss'
                i=0 #color
                num=153
                txt='CG'
            #while(os.path.exists(path+"/rem_iteration"+str(num)+"/CG_hist.txt")):
            #    num+=1
            #num-=1
            #print(num)
                
            x, temp = set_average([num],path)
            for j in range(len(probs_AA)):
                plt.figure(j)
                plt.plot(x[j],np.copy(temp[j]),'-',color=colors[i],linewidth=1,label=txt)
for j in range(len(probs_AA)):
    plt.figure(j)
    plt.xlabel("Distance ($\mathrm{\AA}$)",fontsize=18)
    plt.ylabel("Probability",fontsize=18,labelpad=10)
    plt.legend()
    plt.figure(j).subplots_adjust(left=0.16)
    plt.figure(j).subplots_adjust(bottom=0.15)
    plt.xlim(0,20)
    plt.ylim(0.0,0.0125)
    plt.savefig("site"+str(j)+".png")
    #plt.tight_layout()
    plt.close(j)
    #set_n.append(temp)
         #x, set2_m = set_average([num])
         #x, set3_m = set_average([3])
         #x, set4_m = set_average([4])
            

