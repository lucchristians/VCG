import numpy as np 
from mscg import *

import sys

#find best of 5 virtual sites
def best_fit(vs_can,atoms_,region): #TODO: update this to ensure no repeats if ncessessary
    atoms=np.array(atoms_)
    avg_dist=[]
    min_dist=9999
    min_dist_ind=-1
    for i in range(len(vs_can)):
        
        #print(np.shape(vs_can[i]))
        #print(region)
        #print(np.shape(atroms[region]))
        #print(np.shape(np.linalg.norm(vs_can[i]-pos_array[region],axis=1)))
        #print(np.shape(np.mean(np.linalg.norm(vs_can[i]-pos_array[region],axis=1))))
        avg_dist.append(np.mean(np.linalg.norm(vs_can[i]-atoms[region,:],axis=1)))
        #print(np.shape(avg_dist))
        if(avg_dist[i]<min_dist):
            min_dist=avg_dist[i]
            min_dist_ind=i
    return vs_can[min_dist_ind]

def best_fit_mod(vs_can,atoms_,region,used):
    atoms=np.array(atoms_)
    avg_dist=[]
    min_dist=9999
    min_dist_ind=-1
    for i in range(len(vs_can)):
        done_before=False
        for j in used:
             if(i==j):
                 done_before=True
                 break
        if(done_before):
            continue         
        #print(np.shape(vs_can[i]))
        #print(region)
        #print(np.shape(atroms[region]))
        #print(np.shape(np.linalg.norm(vs_can[i]-pos_array[region],axis=1)))
        #print(np.shape(np.mean(np.linalg.norm(vs_can[i]-pos_array[region],axis=1))))
        avg_dist.append(np.mean(np.linalg.norm(vs_can[i]-atoms[region,:],axis=1)))
        #print(np.shape(avg_dist))
        if(avg_dist[-1]<min_dist):
            min_dist=avg_dist[-1]
            min_dist_ind=i
    used.append(min_dist_ind)
    return used

input_file=sys.argv[1]
output_file=sys.argv[2]

traj=Trajectory(input_file, fmt='lammpstrj')
trajout=Trajectory(output_file, fmt='lammpstrj', mode='w')
while(traj.read_frame()): #TODO: make a stack in python 
    trajout.natoms=int(traj.natoms) + 70
    trajout.box=traj.box
    #define type list
    #trajout.t=np.arange(0,len(traj.t)+60))%(54+6)+1
    #add in stacked vs sites (on bottom chains)
    type_list=[]
    pos_array=[]
    stacked_vs_for_top=[]
    site=[52,53]
    for j in range(len(site)):
        stacked_vs_for_top.append([])
    for i in range(5): #TODO: copy this for top fib
        #trajout.natoms=int(traj.natoms/10 + 1)
        #trajout.box=traj.box
        type_list.append(np.arange(1,62))
        #CG sites
        for j in range(54):
            pos_array.append(traj.x[(i*54+j),:])
        #intraCC vs
        vs_loc=np.array([24,31,38,45,52])-1
        for vs_ in vs_loc:
            pos_array.append(traj.x[(((i+1)%5)*54+vs_),:]) 
        #interCC vs
         #TODO: reverify that this vs is reasonable before proceeding
        
        for jj, j in enumerate(site):
            temp=[]
            for k in range(5):
                temp.append(traj.x[((k+5)*54+j),:]) 
            temp=np.array(temp)
            region=np.array(np.arange(2,10,dtype=int))
        
            pos_array.append(best_fit(temp,pos_array[-59:],region))
            tempe=np.copy(pos_array[-1])
            tempe[0]-=66.0
            stacked_vs_for_top[jj].append(tempe)
    used=[]
    for j in reversed(range(len(site))):
        used.append([])
        
    for i in range(5): #TODO: copy this for top fib
        #trajout.natoms=int(traj.natoms/10 + 1)
        #trajout.box=traj.box
        type_list.append(np.arange(1,62))
        #CG sites
        for j in range(54):
            pos_array.append(traj.x[((i+5)*54+j),:])
        #intraCC vs
        vs_loc=np.array([24,31,38,45,52])-1
        for vs_ in vs_loc:
            pos_array.append(traj.x[((5+(i+1)%5)*54+vs_),:]) 
        #interCC vs
        region=np.array(np.arange(2,10,dtype=int))
        for j in range(len(site)): #TODO: ensure this reflects the order of sites in the lower CC
            used[j] = best_fit_mod(stacked_vs_for_top[j],pos_array[-59:],region,used[j])
        
            #j*=-1
            pos_array.append(stacked_vs_for_top[j][used[j][-1]])
        
    
    pos_array=np.array(pos_array)
    type_list=np.array(type_list)
    
    #trajout.x=traj.x
    trajout.x=pos_array
    trajout.t=type_list
    trajout.write_frame()
    exit()
    
        
        



