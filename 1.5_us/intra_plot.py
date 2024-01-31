import numpy as np
import matplotlib.pyplot as plt
colors=['#000000','#E83274','#4FA3EC']
#colors=['#332288','#117733','#44AA99','#88CCEE','#DDCC77','#CC6677','#AA4499','#882255','#A790CB','#000000']
#load data
aa = np.loadtxt("/mnt/data0/us_results/intra_m2_us/Qwt_us/min_path.dat")
#convert units
aa[:,0]=aa[:,0]*10 # A
aa[:,3]/=4.184 #kJ to kcal

#novs
novs = np.loadtxt("/mnt/data0/cg_pmf_test/novs/intra/all_dat.dat")
#vs
vs = np.loadtxt("/mnt/data0/cg_pmf_test/repl/intra/all_dat.dat")

#plot figures
ratio=3/4
inch=5*4/3
size=[inch,inch*ratio]
plt.figure(figsize=size)
plt.style.use("~/fig_format.mplstyle") #basic formatting
#plt.plot(aa[::5,0],aa[::5,3]-aa[-1,3],'--bo',color=colors[0],label="aa")
#plt.plot(vs[:,0]-4,vs[:,1]-vs[0,1],'--bo',color=colors[1],label="vcg")
#plt.plot(novs[:,0]-4,novs[:,1]-novs[0,1],'--bo',color=colors[2],label="cg")
plt.errorbar(aa[::5,0],aa[::5,3]-aa[-1,3],fmt='--o',color=colors[0],label="AA",markersize=4)
plt.errorbar(novs[:,0]-4,novs[:,1]-novs[0,1],yerr=novs[:,2],fmt='--o',color=colors[1],label="CG",markersize=4)
plt.errorbar(vs[:,0]-4,vs[:,1]-vs[0,1],yerr=vs[:,2],fmt='--o',color=colors[2],label="VCG+HP",markersize=4)
plt.legend()
plt.xlim(4.84,25)
plt.ylim(-30.40,5)
plt.xlabel("Distance ($\mathrm{\AA}$)",labelpad=10)
plt.ylabel("PMF (kcal/mol)",labelpad=10)
plt.tight_layout()
plt.savefig("intra_pmf_fig.png")
print("done")

