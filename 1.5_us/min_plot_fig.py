
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.use('Agg')
data=np.loadtxt("min_path.dat")

plt.figure()
plt.plot(np.mean(data[:,:2],axis=1),data[:,3])
plt.savefig('min_path.png')
