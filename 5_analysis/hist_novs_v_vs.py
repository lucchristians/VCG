import numpy as np
from mscg import *
import matplotlib.pyplot as plt
import os
import sys
 #TODO: create alt version that doesn't do the calculations (use CG_hist.txt)
colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']


max=2.5
L=255.0

path1=sys.argv[1]
path=path1
path2=sys.argv[2]
consts=np.genfromtxt(path1+'/reference_AA/const_AA.txt',ndmin=2)

'''
txt=np.loadtxt('pairhist.txt')
#txt=np.loadtxt('vs_loc')

if(len(np.shape(txt))==1):
	txt = txt.reshape(1,2)
'''
probs_AA = np.loadtxt(path1+"/reference_AA/AA_hist.txt")
#print(np.shape(probs_AA))
if(len(np.shape(probs_AA))==1):
	bins=len(probs_AA)
	probs_AA = probs_AA.reshape(1,bins)
else:
	bins=len(probs_AA[0])

def set_average(path,int_array):
    setn=[]
    x=np.loadtxt(path+"/rem_iteration"+str(int_array[0])+"/CG_ranges.txt")
    for i in int_array:
        #os.chdir(path+"/rem_iteration"+str(i))
        #traj=Trajectory(path+"/rem_iteration"+str(i)+"/sub1.lammpstrj", fmt='lammpstrj') #TODO: name is placeholder
        #x, temp=AA_2_CG(traj,probs_AA,txt,consts)
        temp=np.loadtxt(path+"/rem_iteration"+str(i)+"/CG_hist.txt")
        setn.append(np.copy(temp))
    setn=np.array(setn)

    setn_m=np.mean(setn,axis=0)
    return x, setn_m

'''    
def set_average_split(int_array):
    setn=[]
    setm=[]
    for i in int_array:
        os.chdir(path+"/rem_iteration"+str(i))
        traj=Trajectory("sub1.lammpstrj", fmt='lammpstrj') #TODO: name is placeholder
        x, temp,xe, tempe=AA_2_CG_split(traj,probs_AA,txt,consts)
        setn.append(np.copy(temp))
        setm.append(np.copy(tempe))
    setn=np.array(setn)
    setm=np.array(setm)
    setn_m=np.mean(setn,axis=0)
    setm_m=np.mean(setm,axis=0)
    return x, setn_m, setm_m
'''    
x, set1_m = set_average(path1,[60])
x, set2_m = set_average(path2,[60])
#x, set3_m = set_average([3])
#x, set4_m = set_average([4])
print(np.shape(set1_m))
#print(np.shape(set1_n))

os.chdir(path+"/")
save=True
for i in range(5): #TODO: set to be based on number of vs sites
    name="REM_Distribution_vs%d"%(i)
    fig = plt.figure()
    plt.style.use("~/fig_format.mplstyle") #basic formatting
    plt.plot(x[i],probs_AA[i],'-',color='#CCCCCC',label="AA") #fill
    ax = plt.axes()
    d=np.zeros(len(probs_AA[i]))

    ax.fill_between(x[i],probs_AA[i], where=probs_AA[i]>=d, interpolate=True, color='#CCCCCC')
    
    plt.plot(x[i],set1_m[i],'-',color=colors[i],linewidth=1,label="VCG") ## color blind friendly
    plt.plot(x[i],set2_m[i],'--',color=colors[i],linewidth=1,label="CG")
    #plt.plot(x[i],set3_m[i],'-',color='#EE0000',linewidth=1,label="iter%d"%(7))
    #plt.plot(x[i],set4_m[i],'-',color='#00AAAA',linewidth=1,label="iter%d"%(8)) ## color blind friendly
    #plt.plot(x[i],set1_n[i],'--',color='#00AAFF',linewidth=1,label="iter%d"%(5)) ## color blind friendly
    #plt.plot(x[i],set2_n[i],'--',color='#0000DD',linewidth=1,label="iter%d"%(6))
    #plt.plot(x[i],set3_n[i],'--',color='#EE0000',linewidth=1,label="iter%d"%(7))
    #plt.plot(x[i],set4_n[i],'--',color='#00AAAA',linewidth=1,label="iter%d"%(8))
    #fig.subplots_adjust(left=0.16)
    #fig.subplots_adjust(bottom=0.15)
    plt.xlabel("Distance ($\mathrm{\AA}$)")
    plt.ylabel("Probability")
    #plt.xlim(2.5,9.0)
    plt.legend()
    if(save):
        plt.savefig("%s_equil.png"%(name))
        
    else:
        plt.show()
    plt.close()

