import numpy as np
import matplotlib.pyplot as plt
import sys

vcg=np.loadtxt("sym5_VCG.txt")
cg=np.loadtxt("sym5_CG.txt")
colors=['#000000','#E83274','#4FA3EC']
plt.style.use("~/fig_format.mplstyle") #basic formatting
x=np.arange(18,49)+1
plt.figure()
plt.plot(x,cg[0,:],color=colors[1],label='CG')
plt.fill_between(x,cg[0,:]-cg[1,:],cg[0,:]+cg[1,:],color=colors[1],alpha=0.3)
plt.plot(x,vcg[0,:],color=colors[2],label='VCG')
plt.fill_between(x,vcg[0,:]-vcg[1,:],vcg[0,:]+vcg[1,:],color=colors[2],alpha=0.3)
plt.plot(x,np.mean([cg[0,:],vcg[0,:]],axis=0),color=colors[0],label='mean')
plt.legend()
plt.xlabel("residue")
plt.ylabel("sym score")
plt.savefig('dist_score_comp.png')
print(np.mean(cg[0,:]))
print(np.mean(vcg[0,:]))
print(np.mean([cg[0,:],vcg[0,:]]))
