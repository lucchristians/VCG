import numpy as np
import sys
import matplotlib.pyplot as plt
from mscg import *

def find_potential(name, pot_col_id, block_size):
    distances=[]
    potentials=[] #all runs combined
    skip_first_step=True # if you did minimize
    gather_data=False
    potential=[] #potential per run
    #skip_every_other = False
    skip_every_other = True
    with open(f'{name}/log.pull') as f:
        for line in f:
            data=line.split()
            if(len(data)==0):
                continue
            if(data[0]=='Step'):
                if(skip_first_step):
                    skip_first_step=False
                    continue
                gather_data=True
                #re-initalize frame
                potential=[]
                continue
            
            if(gather_data):
                if(data[0]=='WARNING:') :
                    continue
                elif(data[0]=='Loop'): #process data
                    potential=np.array(potential)

                    if (skip_every_other) :
                        skip_every_other = False
                    else :

                        # perform block analysis
                        ndata = len(potential)
                        #print(ndata)
                        num_blocks = ndata // block_size
                        block_averages = np.zeros(num_blocks)
                        for i in range(num_blocks) :
                            block = potential[i*block_size : (i+1)*block_size]
                            block_averages[i] = np.mean(block)

                        potentials.append(np.copy(block_averages)) # ignore movement time
                        skip_every_other = True
                        
                    gather_data=False
                else :
                    #print(data)
                    potential.append( float(data[pot_col_id])+float(data[pot_col_id-1]) )
    
    potentials=np.array(potentials) # (windows, ndata)
    print(potentials.shape)

    pot_mean = np.mean(potentials, axis=1)
    #pot_std = np.std(potentials, axis=1)
    pot_sem = np.std(potentials, axis=1) / np.sqrt( num_blocks )
    
    #pot_data = np.column_stack( [pot_mean, pot_std] )
    pot_data = np.column_stack( [pot_mean, pot_sem] )
    
    np.savetxt(f'{name}/potentials.dat', potentials)
    return pot_data

def find_dists(name, slice_id) :
    traj = Trajectory(f"{name}/sub1.lammpstrj", fmt='lammpstrj')
    pos_array = []
    cell_lens = traj.box

    while(traj.read_frame()):
        atomtype=traj.t
        pos=traj.x
        pos_array.append(np.copy(pos))
    n_frames=len(pos_array)
    pos_array=np.array(pos_array)

    # first six mols are m1, use COM of types 4-24
    # second six mols are m2

    mask = np.logical_and(atomtype >= 18, atomtype <= 36)
    subset_pos = pos_array[:, mask, :]
    m1_pos = subset_pos[:, :slice_id, :]
    m2_pos = subset_pos[:, slice_id:, :]

    m1_center = np.mean( m1_pos, axis=1 )
    m2_center = np.mean( m2_pos, axis=1 )
    dists = np.linalg.norm( m2_center - m1_center, axis=1 )
    np.savetxt(f'{name}/dists.dat', dists)
    print(dists.shape)
    return dists

#pots=np.zeros([30,0])

#names=[ 'intra_hex', 'inter_hex' ]
#names=['novs/intra','novs/inter','excl/intra','excl/inter','repl/intra','repl/inter']
names=['repl/intra']

# since we have 21 sites counted per monomer, use a factor of 21
slices=[ 19*1 ] 
traj_dt = 1000
log_dt = 50
block_size = 200
for ii, name in enumerate(names) :
    
    #pots = find_potential(name, 4) #4 is the column id with poteng, 10 is the column id with epair
    pots = find_potential(name, 8, block_size) #4 is the column id with poteng, 10 is the column id with epair
    dists = find_dists(name, slices[ii])

    #dists = dists[1:][::10]
    dists = dists[1:]
    #n_dists = len(dists) // 10
    n_dists = len(pots)
    dist_bins = len(dists) // n_dists
    avg_dists = np.zeros( n_dists )
    for jj in range(n_dists) :
        avg_dists[jj] = np.mean( dists[ jj*dist_bins : (jj+1)*dist_bins ] )
    all_dat = np.column_stack( [avg_dists, pots] )
    np.savetxt(f'{name}/all_dat.dat', all_dat)
    

