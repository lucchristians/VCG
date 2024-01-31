import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt("loss.txt",skiprows=1)

plt.figure()
plt.plot(data[:,0],data[:,-1])
plt.savefig("loss.png")
