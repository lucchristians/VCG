import numpy as np
import sys
import matplotlib.pyplot as plt

# process the log file, get the spring forces, convert to displacements
# spring force is f = -k*(r-r0) with e = 0.5*k*(r-r0)^2
# force is negative
def find_disps(kconst, r0, pot_col_id,coord_id):
    distances=[]
    potentials=[] #all runs combined
    skip_first_step=False # if you did minimize
    gather_data=False
    potential=[] #potential per run
    skip_every_other = True
    counter = 1
    with open(f'log.pull') as f:
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
                counter = 1
                continue
            
            if(gather_data):
                if(data[0]=='WARNING:') :
                    continue
                elif(data[0]=='Loop'): #process data
                    potential=np.array(potential)

                    if (skip_every_other) :
                        skip_every_other = False
                    else :
                        potentials.append(potential) 
                        skip_every_other = True
                        
                    gather_data=False
                else :
                    # convert to displacement in A
                    #dx = float(data[pot_col_id])/kconst + r0
                    #dy = float(data[pot_col_id+1])/kconst
                    #dz = float(data[pot_col_id+2])/kconst
                    dcoord=[0,0,0]
                    for i in range(3):
                        #if(data[pot_col_id+i]=='N
                        dcoord[i] = float(data[pot_col_id+i])/kconst + r0
                    
                    #rr = ( dx*dx + dy*dy + dz*dz )**0.5
                    #potential.append( [counter, rr])
                    #potential.append( [counter, dx])
                    potential.append( [counter, dcoord[coord_id]])
                    
                    counter += 1
    
    potentials=np.array(potentials) # (windows, ndata)
    print(potentials.shape)
    # only keep final
    np.savetxt('CV.dat', potentials[-1])
    return potentials

f = open("spring.settings", "r")
flines = f.readlines()
coords=sys.argv[1]
r0_pos={0:7,'x':7,1:8,'y':8,2:9,'z':9}

r0 = float(flines[1].split()[r0_pos[coords]])
kconst = float(flines[1].split()[6])

dict_coords={0:0,'x':0,1:1,'y':1,2:2,'z':2}
coord_ind=dict_coords[coords]
# assumes restraint is [rx, 0, 0]
disps = find_disps(kconst, r0, 12,coord_ind) # 12 is the column index for the x component of the force vector, use all 3 components

    

