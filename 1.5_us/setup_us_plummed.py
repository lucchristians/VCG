import numpy as np
import os
import sys

path=sys.argv[1]
window=int(sys.argv[2])

f = open(path+"us_test.dat",'w') #TODO: change to us.dat once working

#masses for different atoms
atom_masses={'N':14,'O':16,'H':1,'C':12,'S':32}

def read_ndx(file_name):
    group_dict={}
    groupname=''
    with open(file_name) as f:
        for line in f:
            if(len(line)==0):
                continue
            if(line[0]=='['):
                if(groupname!=''):
                    index=np.hstack(index)
                    group_dict[groupname]=index
                groupname=line[line.index('[')+2:line.index(']')-1]
                index=[]
            else:
                index.append(np.array(line.split()).astype(int))
        index=np.hstack(index)
        group_dict[groupname]=index
    return group_dict
                
ndx = read_ndx(path+'../sys.ndx')
#print(ndx)
#exit()
#comp={'x':0,'y':1,'z':2}
comp={0:'x',1:'y',2:'z'}

# interpret data from gro file
with open(path+"d_%d.gro"%(window)) as gro:
    temp_chain=[]
    chain=[]
    pos=[]
    masses=[]
    first_line=True
    second_line=True #TODO: compute mass for each atom
    res_type='Null'
    prev_id=-1
    for line in gro:
        
        if(first_line):
            first_line=False
            continue
        elif(second_line):
            second_line=False
            continue
        
        resid=int(line[:5])
        res_type=line[5:8]
        if(res_type=='SOL'):
            break
        index=int(line[15:20])
        atom=line[11:16].strip()[0]
        x=float(line[20:28])
        y=float(line[28:36])
        z=float(line[36:44])
        #print(resid,res_type,index,x,y,z)
        
        if(resid<prev_id):
            chain.append(temp_chain)
            temp_chain=[]
        temp_chain.append([x,y,z])
        pos.append([x,y,z])
        masses.append(atom_masses[atom])
        prev_id=resid
        #exit()
    chain.append(temp_chain)
    pos=np.array(pos)
    masses=np.array(masses)
    chain=np.array(chain)
    #print(chain[0])
    #exit()    
# shape chains,atoms,coords
#print(np.shape(chain))
n_chains=np.shape(chain)[0] 
chain_len=np.shape(chain)[1]
# define groups for restrants (for variable number of groups)
n_args = int(sys.argv[3]) #must be even
if(n_args%3!=0):
     print("3*n number of arguments need to be specified")
     exit()

#chain_indexs=[]
#groups=[]
group_A=[]
group_B=[]
dists=[]
dists_c=[]
kappa=[]
direction=[]
for i in range(n_args):
    if(i%3==0):
        #TODO: define based on COM of new groups
        '''
        dists.append(np.linalg.norm(np.mean(chain[:5,int(sys.argv[4+i]),:]-chain[5:10,int(sys.argv[5+i]),:],axis=0),axis=0))
        dists_c.append(-1*np.mean(chain[:5,int(sys.argv[4+i]),:]-chain[5:10,int(sys.argv[5+i]),:],axis=0)[comp[sys.argv[6+i]]])
        '''
        #dists.append(np.linalg.norm((np.sum(np.expand_dims(masses[ndx[sys.argv[4+i]]],axis=1)*pos[ndx[sys.argv[4+i]],:],axis=0)/np.sum(masses[ndx[sys.argv[4+i]]]))-(np.sum(np.expand_dims(masses[ndx[sys.argv[5+i]]],axis=1)*pos[ndx[sys.argv[5+i]],:],axis=0)/np.sum(masses[ndx[sys.argv[5+i]]])),axis=0)) #TODO: determine center of mass of this group as needed
        for j in range(3):
            dists_c.append('{dist%s%d}'%(comp[j],int(i/3)))
            if(j==1):
                kappa.append(4500.0)
            else:
                kappa.append(2000.0)
        #dists_c.append(-1*(np.sum(np.expand_dims(masses[ndx[sys.argv[4+i]]],axis=1)*pos[ndx[sys.argv[4+i]],:],axis=0)/np.sum(masses[ndx[sys.argv[4+i]]])-np.sum(np.expand_dims(masses[ndx[sys.argv[5+i]]],axis=1)*pos[ndx[sys.argv[5+i]],:],axis=0)/np.sum(masses[ndx[sys.argv[5+i]]]))[comp[sys.argv[6+i]]])
        #dists_c.append(0.0) #TODO: change this to the above formula for the nonstacked sims
        #print(np.mean(chain[:5,int(sys.argv[4+i]),:]-chain[5:10,int(sys.argv[5+i]),:],axis=0))
        #print(dists[-1])
        #exit()
         #TODO: will this be component based
        direction.append(sys.argv[6+i])
        #print("%7.2f"%(np.max(force[:,1])/2))
        #chain_indexs.append(int(sys.argv[5+i]))
        '''
        group1=[]
        group2=[]
        for j in range(int(n_chains/2)): #TODO: redefine how these groups are defined
            group1.append(int(chain_len*j+int(sys.argv[4+i])))
            group2.append(int(chain_len*(int(n_chains/2)+j)+int(sys.argv[5+i])))
        groups.append([group1,group2])
        '''
        group_A.append(sys.argv[4+i])
        group_B.append(sys.argv[5+i])
#chain_indexs=np.array(chain_indexs)
#groups=np.array(groups)

#dists=np.array(dists)
print(np.shape(direction))
#print(np.shape(dists))
#print(np.shape(groups))
#exit()
# determine Kappa and distance constants for each restraint

# generate file
# groups gen

for i in range(int(n_args/3)):
    f.write("c%d: COM ATOMS="%(2*i+1))
    n=len(ndx[group_A[i]])
    for j in range(len(ndx[group_A[i]])):
        f.write("%d"%(ndx[group_A[i]][j]))
        if(j!=n-1):
            f.write(',')
    '''
    for j in range(n):
        f.write("%d"%(groups[i,0,j]))
        if(j!=n-1):
            f.write(',')
    '''
    f.write('\n')
    f.write("c%d: COM ATOMS="%(2*i+2))
    n=len(ndx[group_B[i]])
    for j in range(len(ndx[group_B[i]])):
        f.write("%d"%(ndx[group_B[i]][j]))
        if(j!=n-1):
            f.write(',')
    '''
    for j in range(n):
        f.write("%d"%(groups[i,1,j]))
        if(j!=n-1):
            f.write(',')
    '''
    f.write('\n')
#dist gen
for i in range(int(n_args/3)): #TODO: set up based on components
    #f.write("dist%d: DISTANCE ATOMS=c%d,c%d\n"%(i,2*i+1,2*i+2))
    f.write("dist%d_c: DISTANCE ATOMS=c%d,c%d COMPONENTS\n"%(i,2*i+1,2*i+2))

#res gen
for i in range(int(n_args/3)):
    #f.write("res%d: RESTRAINT ARG=dist%d AT=%-6.3f KAPPA=%-7.2f\n"%(i,i,dists[i],kappa[i]))
    for j in range(3):
        f.write("res%d_c%s: RESTRAINT ARG=dist%d_c.%s AT=%s KAPPA=%-7.2f\n"%(i,comp[j],i,comp[j],dists_c[3*i+j],kappa[3*i+j]))
f.write("PRINT ARG=")
for i in range(int(n_args/3)):
    #f.write("dist%d"%(i))
    #if(i!=int(n_args/3)-1):
    #f.write(',')
    for j in range(3):
        f.write("dist%d_c.%s"%(i,comp[j]))
        #if(i!=int(n_args/3)-1):
        f.write(',')
for i in range(int(n_args/3)):
    for j in range(3):
        f.write("res%d_c%s.bias"%(i,comp[j]))
        if(i!=int(n_args/3)-1 or j!=2):
            f.write(',')
f.write(" STRIDE=5000 FILE=colvar")
    
exit()
#get c1 and c2 
c1=''
c2=''
#read data
#with open("../sys.ndx") as data:
#    ATE_reached=False
#    FTJ_reached=False
#    for line in data:
#        #TODO: use thsi as a range not print groups
#        if(line=='[ ATE ]\n'):
#            ATE_reached=True
#            continue
#        if(line=='[ FTJ ]\n'):
#            FTJ_reached=True
#            continue
#
#       if(FTJ_reached):
#            temp=line.split()
#            for index in temp:
#                c2+=str(index)+","
#            continue
#        elif(ATE_reached):
#            temp=line.split()
#            for index in temp:
#                c1+=str(index)+","
#            continue
#print(c1[:-1])
#print(c2[:-1])

#determine distance
dist=np.loadtxt('../summary_distances.dat')
x0=dist[0,1]
dx=0.3
path = os.getcwd()
index = int(path.split('/')[-1][-1]) #TODO may not be needed anymore
print("%6.3f"%(x0+index*dx))

#determine kappa
force=np.loadtxt("../pull_pullf.xvg",skiprows=17)
print("%7.2f"%(np.max(force[:,1])/2))
