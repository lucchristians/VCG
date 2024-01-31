import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as st
colors=['#000000','#E83274','#4FA3EC']
paths=['gauss/logs/','vsite/logs/']
#colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']
#load in data
path1=paths[0]#sys.argv[1]
path2=paths[1]#sys.argv[2]
loss1=[]
for i in range(4):
    loss1.append(np.loadtxt(path1+r'loss_mse_%d.txt'%(i+1)))
loss1=np.array(loss1)
print(loss1.shape)
mean_loss1=np.mean(np.mean(loss1[:,:,1:-3],axis=2),axis=0)
stde_loss1=np.std(np.mean(loss1[:,:,1:-3],axis=2),axis=0)
print(mean_loss1.shape)
loss2=[]
for i in range(4):
    loss2.append(np.loadtxt(path2+r'loss_mse_%d.txt'%(i+1)))
loss2=np.array(loss2)
mean_loss2=np.mean(np.mean(loss2[:,:,1:-3],axis=2),axis=0)
stde_loss2=np.std(np.mean(loss2[:,:,1:-3],axis=2),axis=0)

#ttest

#all losses
all_ttest= st.ttest_ind(np.mean(loss1[:,:300,1:-3],axis=2),np.mean(loss2[:,:300,1:-3],axis=2))
#min losses
loss1_min=np.min(np.mean(loss1[:,:,1:-3],axis=2),axis=1)
loss2_min=np.min(np.mean(loss2[:,:,1:-3],axis=2),axis=1)
print(loss1_min,loss2_min)
min_ttest= st.ttest_ind(loss1_min,loss2_min)
#print("all loss t-test")
#print(all_ttest)
vmin_ind=243
cmin_ind=153
min_vloss=mean_loss2[vmin_ind]
min_vloss_std=stde_loss2[vmin_ind]
min_closs=mean_loss1[cmin_ind]
min_closs_std=stde_loss1[cmin_ind]


print("min loss t-testt")
print(min_ttest)
print('VCG')
print(min_vloss,min_vloss_std)
print('CG')
print(min_closs,min_closs_std)

exit()

#plt.style.use("~/fig_format.mplstyle") #basic formatting

ratio=3/4
inch=5*4/3
size=[inch,inch*ratio]

plt.figure(figsize=size)
plt.style.use("~/fig_format.mplstyle") #basic formatting
plt.plot(loss1[0,:,0],mean_loss1,color=colors[1],label='CG')
plt.fill_between(loss1[0,:,0],mean_loss1-stde_loss1,mean_loss1+stde_loss1,color=colors[1],alpha=0.3)
plt.plot(loss2[0,:,0],mean_loss2,color=colors[2],label='VCG')
plt.fill_between(loss2[0,:,0],mean_loss2-stde_loss2,mean_loss2+stde_loss2,color=colors[2],alpha=0.3)
#np.savetxt('loss_CG.txt',np.mean(loss1[:,1:-3],axis=1))
#np.savetxt('loss_VCG.txt',np.mean(loss2[:,1:-3],axis=1))
plt.legend()
plt.xlim(0,300)
plt.ylim(0.0002,0.0003)
plt.xlabel("Iteration",labelpad=10)
plt.ylabel("MSE Loss",labelpad=10)
plt.tight_layout()
plt.savefig("REM_loss.png")
