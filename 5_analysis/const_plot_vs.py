import numpy as np
import matplotlib.pyplot as plt
import sys

path=sys.argv[1]
if(path[-1]!='/'):
    path+='/'
#system=path.split('_')[2]
system='novs' 
#universal parameters
colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']
fmt=['-','--','-.']

xlabel="Iteration"
if(system=='vs'):
    ylabel=["K","A","B"]
else:
    ylabel=["A","B"]
#fnt_size=17
#lgnd_size=8
data=np.loadtxt(path+"constants.txt")
points_per_const=int(np.shape(data)[1]/len(ylabel))
for i in range(len(ylabel)):
	title="%s_const"%(ylabel[i])
	
	plt.figure()
	plt.style.use("~/fig_format.mplstyle") #basic formatting
	for j in range(points_per_const):
		plt.plot(data[:,0],data[:,1+i+j*len(ylabel)],'-',color=colors[j%len(colors)],label=ylabel[i]+str(j))
	#plt.xlim(xlim)
	#plt.ylim(ylim)
	plt.tick_params(axis="x",direction="in")
	plt.tick_params(axis="y",direction="in")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel[i])
	plt.legend()
	plt.tight_layout()
	#plt.title(title,fontsize=fnt_size+3)
	plt.savefig(path+title+"_plot",dpi=300)

