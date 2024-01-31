import numpy as np
import matplotlib.pyplot as plt

colors=['#000000','#E83274','#4FA3EC']

rmsd=[]

for i in range(4):
    temp=np.loadtxt('run%d/rmsd_CA.xvg'%(i),skiprows=18)
    x=np.copy(temp[:,0])
    rmsd.append(temp[:,1])
'''
for i in range(4):
    temp=np.loadtxt('run%d/rmsd_top.xvg'%(i),skiprows=18)
    x=np.copy(temp[:,0])
    rmsd.append(temp[:,1])
    temp=np.loadtxt('run%d/rmsd_bot.xvg'%(i),skiprows=18)
    x=np.copy(temp[:,0])
    rmsd.append(temp[:,1])
'''
rmsd=np.array(rmsd)
length= 660 #ns
x=np.linspace(0,length,num=len(x))

mean_rmsd=np.mean(rmsd,axis=0)
stde_rmsd=np.std(rmsd,axis=0)

ratio=3/4
inch=5*4/3
size=[inch,inch*ratio]

plt.figure(figsize=size)
plt.style.use("~/fig_format.mplstyle") #basic formatting
plt.plot(x,mean_rmsd,color=colors[1],label='CG')
plt.fill_between(x,mean_rmsd-stde_rmsd,mean_rmsd+stde_rmsd,color=colors[1],alpha=0.3)
plt.xlabel("Time (ns)",labelpad=10)
plt.ylabel("RMSD (nm)",labelpad=10)
plt.xlim(0,200)
plt.ylim(0,1.4)
#plt.show()
plt.savefig('rmsd_mean.png')
