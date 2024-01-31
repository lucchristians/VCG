import numpy as np
import sys
import matplotlib.pyplot as plt

array_len=len(sys.argv)-1

data=[]
print("path | loss | dloss | min_index")

for i in range(array_len):
    #print(sys.argv[i+1])
    temp=np.loadtxt(sys.argv[i+1],skiprows=1)
    plt.figure()
    path=sys.argv[i+1][:sys.argv[i+1].index('/')+1]
    plt.plot(temp[:,0],temp[:,-1])
    plt.savefig(path+"loss.png")
    plt.close()
    loss=min(temp[:,-1]) #last column; last row
    loss_ind=np.argmin(temp[:,-1])
    dloss=np.copy(temp[0,-1]-min(temp[:,-1]))
    data.append([path,loss,dloss])

    print(path,loss,dloss,loss_ind)


