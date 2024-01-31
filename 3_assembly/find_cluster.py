import numpy as np
import networkx as nx
import sys as sys
import math as math
from mscg import *

# USAGE:
# analyze_assembly_and_extract_clusters.py trajbasename startindex endindex clusterMolID dtau_per_frame 
# this version assumes trajectory is LAMMPSTRJ format custom (x, y, z)

##################### USER SPECIFIED INFO ##############################

path=sys.argv[1]
trajname=sys.argv[2] # The name of the trajectory (base+index.lammpstrj)
out_ind=int(sys.argv[3])
print("starting ",path,trajname,out_ind)
##################### END USER SPECIFIED INFO ##############################

# we will check
# inter: 3-53 (r=10.76)
# intra: 37-38 (r=5.06)
cont_dist1 = 0.0 #12.0 #distance for contacts (real-real for inter1)
cont_dist2 = 7.0 #distance for contacts (real-real for intra1)

FSTATS=open(path+"stas%d.dat"%(out_ind),"w")
FAST_FIND_BIN = 30.0
DEBUG_FLAG = 0


#ncg = 61 #number of cg sites per mol

FSTATS.write("#time cluster1_size\n")

# write lammps stuff
def write_lammps_header(tfile, time, natoms, dims, fields) :
    tfile.write("ITEM: TIMESTEP\n")
    tfile.write("%d\n" % time)
    tfile.write("ITEM: NUMBER OF ATOMS\n")
    tfile.write("%d\n" % natoms)
    tfile.write("ITEM: BOX BOUNDS pp pp pp\n")
    tfile.write("0.0 %f\n" % dims[0])
    tfile.write("0.0 %f\n" % dims[1])
    tfile.write("0.0 %f\n" % dims[2])
    tfile.write("ITEM: ATOMS ")
    for f in fields :
        tfile.write("%s " % f)
    tfile.write("\n")

def write_lammps_atom(tfile, index, type, x, y, z, vx, vy, vz) :
    tfile.write("%d %d %f %f %f %f %f %f\n" % (index, type, x, y, z, vx, vy, vz))

#returns array ( nframes, natoms, (aid,type,x,y,z) )
def load_atom_traj(trajname) : 
    i=[]
    t=[]
    pos_array = []
    #print(i)
    traj=Trajectory(trajname, fmt='lammpstrj')
    num=-1
    while(traj.read_frame()):
        num+=1
        if(num%10==0):
            #print(num,end='\r')
            atomtype=traj.t
            box=traj.box
            i.append(np.arange(0,len(atomtype),1).astype(int).reshape(-1,1))
            t.append(atomtype.reshape(-1,1))
            pos_array.append(np.copy(traj.x))
    #print("loaded")
    nframes=len(pos_array)
    natoms=len(atomtype)
    i=np.array(i)
    t=np.array(t)
    pos_array=np.array(pos_array)
    it=np.append(i,t,axis=2)
    traj_array=np.append(it,pos_array,axis=2)
    
    return traj_array, box, nframes, natoms

def calc_dist(r1, r2, dims) :
    dr = r2 - r1
    dr = dr - dims * np.rint( dr / dims )
    '''
    dx = abs(r2[0] - r1[0])
    dy = abs(r2[1] - r1[1])
    dz = abs(r2[2] - r1[2])
    if(dx > dims[0]/2.0) :
        dx - dims[0]
    if(dy > dims[1]/2.0) :
        dy - dims[1]
    if(dz > dims[2]/2.0) :
        dz - dims[2]
    '''
    
    dist = np.linalg.norm(dr)

    return dist

##################### COLLECT ALL TRAJECTORY DATA FIRST ######################################
nframes = 0
traj, dims, nframes, natoms = load_atom_traj(path+"%s.lammpstrj" % trajname)
ncg=int(natoms/484)
print(traj.shape)


##################### ANALYZE ALL TRAJECTORY DATA ######################################
# analyze cluster size

outfile = open(path+"intra%d.lammpstrj"%(out_ind),"w")

for i in range(nframes) :

    int1_A_indices = set() #interface1 - sideA
    int1_B_indices = set() #interface1 - sideB
    int2_A_indices = set() #interface2
    int2_B_indices = set() #interface2

    # find all indices of int1 and int2 cg types
    for ff, f in enumerate(traj[i,:,1]) : #gets cg type
        if(f == 3) :
            int1_A_indices.add(ff)
        if(f == 53) :
            int1_B_indices.add(ff)
        if(f == 37) :
            int2_A_indices.add(ff)
        if(f == 38) :
            int2_B_indices.add(ff)
    # two interfaces 3_53 (inter) and 37_38 (intra)
     
    nmol = len(int1_A_indices) #num mols

    # fast neighbor search
    nbinX = int( dims[0] / FAST_FIND_BIN )
    nbinY = int( dims[1] / FAST_FIND_BIN )
    nbinZ = int( dims[2] / FAST_FIND_BIN )
    
    # limits number of interactions to calculate to spped up code
    fast_find_grid = np.zeros((nbinX, nbinY, nbinZ, natoms+1)) #0-index has length of list, rest is list of indices in grid
    for ff, (x,y,z) in enumerate(traj[i,:,2:5]) :
        ffx = int(x/FAST_FIND_BIN)
        ffy = int(y/FAST_FIND_BIN)
        ffz = int(z/FAST_FIND_BIN)
        if(ffx < 0) :
            ffx = ffx + nbinX
        if(ffx >= nbinX) :
            ffx = ffx - nbinX
        if(ffy < 0) :
            ffy = ffy + nbinY
        if(ffy >= nbinY) :
            ffy = ffy - nbinY
        if(ffz < 0) :
            ffz = ffz + nbinZ
        if(ffz >= nbinZ) :
            ffz = ffz - nbinZ
        ffcount = int(fast_find_grid[ffx,ffy,ffz,0])
        fast_find_grid[ffx,ffy,ffz,0] = ffcount + 1
        fast_find_grid[ffx,ffy,ffz,ffcount+1] = ff
    
    # process connectivity graph
    print("Processing graph for frame %d" % i)
    G = nx.Graph() #gag-gag
    for m in range(nbinX) :
        for n in range(nbinY) :
            for l in range(nbinZ) :
                nn1_num = int(fast_find_grid[m,n,l,0])
                for nn1_id in range(nn1_num) :
                    aid1 = int(fast_find_grid[m,n,l,nn1_id+1])
                    #check within +-1 of the current grid
                    for mm in range(3) :
                        for nn in range(3) :
                            for ll in range(3) :
                                mtemp = m - 1 + mm
                                ntemp = n - 1 + nn
                                ltemp = l - 1 + ll
                                if(mtemp < 0) :
                                    mtemp = mtemp + nbinX
                                if(mtemp >= nbinX) :
                                    mtemp = mtemp - nbinX
                                if(ntemp < 0) :
                                    ntemp = ntemp + nbinY
                                if(ntemp >= nbinY) :
                                    ntemp = ntemp - nbinY
                                if(ltemp < 0) :
                                    ltemp = ltemp + nbinZ
                                if(ltemp >= nbinZ) :
                                    ltemp = ltemp - nbinZ
                                nn2_num = int(fast_find_grid[mtemp,ntemp,ltemp,0])
                                for nn2_id in range(nn2_num) :
                                    aid2 = int(fast_find_grid[mtemp,ntemp,ltemp,nn2_id+1])
                                    # now check for contacts
                                    if(aid1 == aid2) : #ignores if from the same molecule
                                        continue
                                    if(aid1 in int1_A_indices and aid2 in int1_B_indices) :
                                        r1 = traj[i, aid1, 2:5]
                                        r2 = traj[i, aid2, 2:5]
                                        m1 = int(math.ceil( float(traj[i, aid1, 0]) / float(ncg) ) )
                                        m2 = int(math.ceil( float(traj[i, aid2, 0]) / float(ncg) ) )
                                        dist = calc_dist(r1, r2, dims)
                                        #print("in first check")
                                        if ( dist < cont_dist1 ) : #checks distance to see if bond exists (inter)
                                            G.add_edges_from([(m1, m2)])
                                            #print("adding a 1-1 edge between: ", m1, m2)
                                    elif(aid1 in int2_A_indices and aid2 in int2_B_indices) :
                                        r1 = traj[i, aid1, 2:5]
                                        r2 = traj[i, aid2, 2:5]
                                        m1 = int(math.ceil( float(traj[i, aid1, 0]) / float(ncg) ) )
                                        m2 = int(math.ceil( float(traj[i, aid2, 0]) / float(ncg) ) )
                                        dist = calc_dist(r1, r2, dims)
                                        if ( dist < cont_dist2 ) : #checks distance to see if bond exists (intra)
                                            G.add_edges_from([(m1, m2)])
                                            #print("adding a 2-2 edge between: ", m1, m2)

    cluster_size = 0.0 #equivalent to nc_gag
    total_size = 0.0
    nc_mol = 0.0
    G_max = set() #will be a set

    try : 
        #G_max = nx.node_connected_component(G, mid)
        Graph_s = [G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=True)]
        
        Graph_lens = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
        cluster_size=Graph_lens[0]
        average_size=np.mean(Graph_lens)
        print(Graph_lens)
        print(average_size)
        print(Graph_lens[0])
        
    except Exception as e:
        print("No cluster found")
        print(e)

    time = i

    # only continue if we have a cluster size > 1, otherwise, skip printing?
    if ( cluster_size > 0 ) :
        num=0 
        for Graph in Graph_s:
            num+=1
            if(num>5):
                break
            for nn in Graph : #TODO: modify this to work for not just the largest graph maybe
                if(nn <= nmol) :
                    nc_mol += 1.0
            # print padded trajectory with clustered mols
            if(num>1):
                continue
            tempatoms = int(nc_mol*ncg)
            write_lammps_header(outfile, time, tempatoms, dims, ["id", "type", "x", "y", "z", "vx", "vy", "vz"])
            for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :
                
                # get molID then determine if its in the GRI_max cluster
                # if so, print it
                tempmid = 0
                if ( aid <= nmol*ncg) :
                    tempmid = int(math.ceil( aid / float(ncg) ))
                if ( tempmid in Graph ) : #TODO: changed from G_max
                    write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, num, 0)
                        

    FSTATS.write("%d %f\n" % (time, nc_mol))        

FSTATS.close()
outfile.close()
