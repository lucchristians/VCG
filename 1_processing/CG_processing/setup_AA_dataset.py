import numpy as np
import matplotlib.pyplot as plt

vs_intra=np.loadtxt("vs_r0_stacked/AA_hist.txt")
#vs_min=np.loadtxt("rcut_stacked/CG_hist_min.txt")[5]
vs_subset=np.loadtxt("rcut_stacked/CG_hist.txt")[:3]
#vs_new_min=[]
vs_new_subset=[]

for i in range(len(vs_intra)):
    #vs_new_min.append(vs_intra[i]/2)
    vs_new_subset.append(vs_intra[i]/2)

#vs_new_min.append(vs_min)
for i in range(len(vs_subset)):
    vs_new_subset.append(vs_subset[i])

#np.savetxt("vs_r0_stacked/AA_hist_min.txt",vs_new_min)
np.savetxt("vs_r0_stacked/AA_hist_subset.txt",vs_new_subset)


