import numpy as np
import matplotlib.pyplot as plt
import sys

paths=['gauss/','vsite/']
#if(path[-1]!='/'):
#    path+='/'
#system=path.split('_')[2]
system=['novs','vs'] 
#universal parameters
ylim=np.array([[[-9,-3],[0.55,1.1]],[[2,12],[0,1.1]]])
#yticks=np.array([[[-8,-6,-4],[0.6,0.8,1.0]],[[4,7,10],[0.2,0.6,1.0]]])
yticks=np.zeros([2,2,5])
for i in range(2):
    for j in range(2):
        p1=ylim[i,j,0]+(ylim[i,j,1]-ylim[i,j,0])*1/4
        p2=ylim[i,j,0]+(ylim[i,j,1]-ylim[i,j,0])*2/4
        p3=ylim[i,j,0]+(ylim[i,j,1]-ylim[i,j,0])*3/4
        
        yticks[i,j,:]=np.array([ylim[i,j,0],p1,p2,p3,ylim[i,j,1]])
#colors=['#000000','#E83274','#4FA3EC']
colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']
#fmt=['-','--','-.']
ratio=2.2
inch=10/2
size=[inch*ratio,inch]

plt.style.use("~/fig_format.mplstyle") #basic formatting
fig, axs = plt.subplots(nrows=2,ncols=2,sharex=True,figsize=size)

plt.subplots_adjust(hspace=0.13,wspace=0.5,right=0.85)
#axs.yaxis.set_major_formatter('{y:9<5.2f}')
for t, sys_ in enumerate(system):
    path=paths[t]
    xlabel="Iteration"
    if(sys_=='vs'):
        ylabel=["K","A (kcal/mol)","B ($\mathrm{\AA}^{-2}$)"]
        
    else:
        ylabel=["H (kcal $\mathrm{\AA}$/mol)",r"$\sigma$ ($\mathrm{\AA}$)"]
    simple_leg='VS'
    #fnt_size=17
    #lgnd_size=8
    data=np.loadtxt(path+"constants.txt")
    points_per_const=int(np.shape(data)[1]/len(ylabel))
    for i in range(len(ylabel)):
        if(sys_=='vs'):
            offset=-1
            if(i==0):
                print("happy")
                continue
        else:
            offset=0
        #title="%s_const"%(ylabel[i])
            
        #plt.figure()
        #plt.style.use("~/fig_format.mplstyle") #basic formatting
        for j in range(points_per_const-2):
            axs[i+offset,t].plot(data[:,0],data[:,1+i+j*len(ylabel)],'-',color=colors[j%len(colors)],label=simple_leg+str(j+1))
        axs[i+offset,t].set_xlim(-0.1,300)
        axs[i+offset,t].set_yticks(yticks[t,i+offset])
        axs[i+offset,t].set_ylim(ylim[t,i+offset])
        #plt.tick_params(axis="x",direction="in")
        #plt.tick_params(axis="y",direction="in")
        if(i+offset==1):
            axs[i+offset,t].set_xlabel(xlabel,labelpad=7)
        #if(t==0):
        axs[i+offset,t].set_ylabel(np.copy(ylabel[i]),labelpad=10)
axs[1,1].legend(loc=(1.05,0.65))
fig.align_ylabels(axs[:, 0])
fig.align_ylabels(axs[:, 1])
#lt.tight_layout()
        #plt.title(title,fontsize=fnt_size+3)
plt.savefig("REM_const.png")

