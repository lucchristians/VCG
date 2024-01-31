import numpy as np
from mscg import *
import matplotlib.pyplot as plt
import os

max=2.5
L=255.0
path = os.getcwd()
def MSE_loss(AA_hist,CG_hist):
    return np.average((AA_hist-CG_hist)**2)

consts=np.genfromtxt(path+'/reference_AA/const_AA.txt',ndmin=2)

path = os.getcwd() 
AA_hist = np.loadtxt(path+"/reference_AA/AA_hist.txt")
#print(np.shape(probs_AA))
if(len(np.shape(AA_hist))==1):
	bins=len(AA_hist)
	AA_hist = AA_hist.reshape(1,bins)
#print(AA_hist)
num=1
iteration=[]
loss=[]
while(os.path.exists(path+'/rem_iteration%d'%(num))):
    temp=0
    CG_hist=np.loadtxt(path+"/rem_iteration"+str(num)+"/CG_hist.txt")
    for i in range(len(AA_hist)):
        temp += MSE_loss(AA_hist[i],CG_hist[i])
    iteration.append(num)
    loss.append(temp)
    num+=1
iteration=np.array(iteration)
loss=np.array(loss)

np.savetxt(path+"/mse_loss.txt",loss)

