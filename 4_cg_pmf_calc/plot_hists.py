import numpy as np
import sys 
import os
import matplotlib.pyplot as plt

var_type=sys.argv[1]

def find_r0(fName,var_):
    var_name='res'+var_[4:].replace('.','')
    with open(fName) as f:
        for line in f:
            if(len(line)==0):
                continue
            if(line[0]=='c' or line[0]=='d' or line[0]=='P'):
                continue
            elif(line[:line.index(':')]==var_name):
                data=line.split()
                num=float(data[3][data[3].index('=')+1:])
                break
    return num

def find_ci(fName,var_):
    with open(fName) as f:
        for line in f:
            data = line.split()
            for i in range(len(data)):
                if(data[i]==var_):
                    c_ind=i-2
                    break
            break
    return c_ind


r0_dists=[]
paths=[]
for i in range(len(sys.argv)-2):
    paths.append(sys.argv[i+2])
    #print(sys.argv[i+2])
    r0_dists.append(find_r0(sys.argv[i+2]+'/us.dat',var_type))
#del temp


plt.style.use("/home/lchristians/fig_format.mplstyle") #basic formatting
size=1.8
h=2.04*size
w=3.15*size
plt.figure(figsize=(w,h))
plt.subplots_adjust(top=0.8)
for i in range(len(paths)):
    #print(paths[i])
    data=np.loadtxt(paths[i]+'/colvar',skiprows=1) #TODO: get restraints and isolate different files
    column_ind=find_ci(paths[i]+'/colvar',var_type)
    hist, x=np.histogram(data[:,column_ind],range=(0.5,3.5),bins=1000,density=True)
    plt.plot([r0_dists[i],r0_dists[i]],[0,max(hist)],'--',color='#797979',linewidth=2)
    plt.plot(np.copy(x[1:]),np.copy(hist),linewidth=2)

    #plt.plot(data[:,0],data[:,column_ind])
plt.xlabel("Distance (nm)")
plt.ylabel("Probability")
plt.xlim(0.5,2.5)
plt.ylim(0,30)
plt.tight_layout()
plt.savefig("hists.png",dpi=300)
