import numpy as np
import matplotlib.pyplot as plt
import sys
colors_t=['#CCCCCC','#E83274','#4FA3EC']
colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']
markers=['o','s','v','D','p','h']
path1=sys.argv[1]
path2=sys.argv[2]
outputName=sys.argv[3]
vs_rep=sys.argv[4]
if(path1[-1]=='/'):
    path1=path1[:-1]
if(path2[-1]=='/'):
    path2=path2[:-1]
def get_sizes(path):
    rep_sizes=[]
    for i in range(4):
        sizes=[]
        with open(path+'_%d/'%(i+1)+'sizes.txt') as f:
            print(i)
            read_data=False
            #sizes=[]
            for line in f:
                if(line.split()[0]=='Processing'):
                    read_data=True
                    continue
                if(read_data):
                    data=line.split(', ')
                    data[0]=data[0][1:]
                    data[-1]=data[-1].strip()[:-1]
                    data=np.array(data).astype(int)
                    s=[]
                    for i in range(8):
                        size=i+1
                        count=len(data[data==size])
                        s.append(count)
                    sizes.append(s)
                    read_data=False
                    continue
            rep_sizes.append(np.copy(np.array(sizes).T))
    rep_sizes=np.array(rep_sizes)
    print(rep_sizes.shape)
    avg_sizes=np.mean(rep_sizes,axis=0)
    std_sizes=np.std(rep_sizes,axis=0)/np.sqrt(rep_sizes.shape[0])
    total=np.sum(avg_sizes,axis=0)
    counts=2*rep_sizes[:,1,:]+3*rep_sizes[:,2,:]+4*rep_sizes[:,3,:]+5*rep_sizes[:,4,:]
    counts_avg=np.mean(counts,axis=0)
    counts_std=np.std(counts,axis=0)
    t=np.arange(0,len(counts_avg))
    print(np.max(std_sizes))
    print(np.max(counts_std))
    return t, avg_sizes, std_sizes, total, counts_avg, counts_std
t1, sizes1, sem1, total1, counts1, csem1 = get_sizes(path1)
t2, sizes2, sem2, total2, counts2, csem2 = get_sizes(path2)

ratio=3/7
inch=5*7/3
size=[inch+1,inch*ratio]
plt.style.use("~/fig_format.mplstyle") #basic formatting
fig = plt.figure(figsize=size)
gs = fig.add_gridspec(1,2,wspace=0.075)
ax = gs.subplots(sharex=True, sharey=True)
#fig, ax = plt.subplots(1,2,figsize=size) #TODO: either create two figures or 2 subplots
#fig.subplots_adjust(right=0.8)
plot_lines=[]
for i in range(len(sizes1)-2):
    #plt.fill_between(t1,sizes1[i]/total1-sem1[i]/total1,sizes1[i]/total1+sem1[i]/total1,color=colors[i],alpha=0.3)
    l1 = ax[0].errorbar(t1[::10],sizes1[i][::10]/total1[::10],yerr=sem1[i][::10]/total1[::10],linestyle='-',label='%d-mer'%(i+1),color=colors[i],markersize=4,marker=markers[i])
    #l1, = plt.errorbar(t1,sizes1[i]/total1,yerr=sem1[i]/total1,linestyle='-',label='%d-mer'%(i+1),color=colors[i],markersize=4,marker=markers[i])
    #plt.fill_between(t2,sizes2[i]/total2-sem2[i]/total2,sizes2[i]/total2+sem2[i]/total2,color=colors[i],alpha=0.3)
    l2 = ax[1].errorbar(t2[::10],sizes2[i][::10]/total2[::10],yerr=sem2[i][::10]/total2[::10],linestyle='-',label='%d-mer'%(i+1),color=colors[i],markersize=4,marker=markers[i])
    plot_lines.append([l1,l2])
ax[1].label_outer()
#legend1 = plt.legend(plot_lines[0], ['CG',vs_rep],loc=(1.00,0.32)) # .79 .82
#fig.gca().add_artist(legend1)
legend2 = plt.gca().legend([l[0] for l in plot_lines], ['%d-mer'%(i+1) for i in range(len(plot_lines))],loc=(1.00,0.535))
fig.gca().add_artist(legend2)
#plt.plot(counts1)
#plt.title(path1[:-1])
#plt.ylim(0,20)
ax[0].set_xlim(0,400)
ax[0].set_ylim(0,1)
ax[1].set_xlim(0,400)
ax[1].set_ylim(0,1)
ax[0].set_xlabel(r"Time ($\tau$)",labelpad=0,fontname='Arial')
ax[1].set_xlabel(r"Time ($\tau$)",labelpad=0,fontname='Arial')
ax[0].set_ylabel(r'$N_{oligomer}/N_{tot}$',labelpad=10,fontname='Arial')
plt.gcf().text(0.462,0.830,r'CG',fontsize=18,fontname='Arial')
plt.gcf().text(0.807,0.830,r'VCG+HP',fontsize=18,fontname='Arial') #TODO: update naming convention

plt.gcf().text(0.482,0.030,r'$x10^{6}$',fontsize=14,fontname='Arial')
plt.gcf().text(0.883,0.030,r'$x10^{6}$',fontsize=14,fontname='Arial')
#plt.xlim(0,30)
#plt.show()
plt.savefig(outputName)


plt.figure(figsize=size)
plt.style.use("~/fig_format.mplstyle") #basic formatting
plt.fill_between(t1,counts1-csem1,counts1+csem1,color=colors_t[1],alpha=0.3)
plt.plot(t1,counts1,'-',label="CG",color=colors_t[1],markersize=4)
plt.fill_between(t2,counts2-csem2,counts2+csem2,color=colors_t[2],alpha=0.3)
plt.plot(t2,counts2,'-',label=vs_rep,color=colors_t[2],markersize=4)
text1=plt.gcf().text(0.869,0.030,'$x10^{6}$',fontsize=14,fontname='Arial')
#plt.title(path1[:-1])
leg1 = plt.legend()
#plt.ylim(0,20)
plt.xlim(0,400)
plt.ylim(0,500)
plt.xlabel(r"Time ($\tau$)",labelpad=7)
plt.ylabel("# of Q in Oligomers",labelpad=10)
plt.savefig("total_"+outputName)
plt.xlim(0,10)
plt.ylim(0,300)
plt.xlabel(None,labelpad=7)
plt.ylabel(None,labelpad=10)
plt.xticks([0.0,5,10],fontsize=25)
plt.yticks([0,150,300],fontsize=25)
leg1.remove()
text1.remove()

plt.savefig("zoomed_"+outputName)

