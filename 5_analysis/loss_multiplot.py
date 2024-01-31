import numpy as np
import matplotlib.pyplot as plt
import sys

pathName=sys.argv[1]
count=int(sys.argv[2])

#name=path[:-1]
plt.figure()
for i in range(count):
    temp = np.loadtxt(pathName+str(i+1)+'/loss.txt')
    plt.plot(temp[:,0],temp[:,-1],label='rep'+str(i))
plt.xlabel("Iterations")
plt.ylabel("L2 Loss")
plt.title(pathName[:-1])
plt.legend()
plt.savefig(pathName[:-1]+'.png')

