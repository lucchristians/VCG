import numpy as np
import matplotlib.pyplot as plt
import sys

colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']

#load in data
path1=sys.argv[1]
path2=sys.argv[2]
loss1=np.loadtxt(path1+'mse_loss.txt')
loss2=np.loadtxt(path2+'mse_loss.txt')

plt.figure()
plt.style.use("~/fig_format.mplstyle") #basic formatting
plt.plot(loss1[:,0],loss1[:,-1],color=colors[0],label='CG')
plt.plot(loss2[:,0],loss2[:,-1],color=colors[1],label='VCG')
plt.legend()
plt.xlim(0,300)
plt.ylim(0.001,0.0045)
plt.xlabel("Iteration")
plt.ylabel("MSE Loss")
plt.tight_layout()
plt.savefig(path2+"intra_novs_v_vs_loss.png")
