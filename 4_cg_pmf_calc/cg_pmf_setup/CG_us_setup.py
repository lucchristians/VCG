import numpy as np
import sys

#initalize group directories
groupTypes={'a':'id','t':'type','m':'molecule'}
dir={'x':0,0:0,'y':1,1:1,'z':2,2:2}
def groupinterrepter(string):
    temp=string[1:].replace(',',' ')

    return "%s %s"%(groupTypes[string[0]],temp)

outputName=sys.argv[1]
direction=dir[sys.argv[2]] # 0:x 1:y 2:z
groupFix=groupinterrepter(sys.argv[3])
n_groups=int(sys.argv[4])
bins=int(sys.argv[5])
groups=[]
r0=[]
rf=[]
k=[]
for i in range(n_groups):
    base=6+i*5
    group1=groupinterrepter(sys.argv[base]) #fmt: a1, m1:5, a-atom, m-molecule, t-type, :- within range (inclusive), ,- different groups 
    group2=groupinterrepter(sys.argv[base+1])
    groups.append([group1,group2])
    r0.append(float(sys.argv[base+2]))
    rf.append(float(sys.argv[base+3]))
    k.append(sys.argv[base+4])

dr=[]
for i in range(len(r0)):
    dr.append((rf[i]-r0[i])/(bins-1))

f=open(outputName,'w')
#define restraints
f.write("group res %s\n"%(groupFix))
f.write("fix vsi res move linear 0 0 0\n")

#interpret & define groups
for i in range(len(groups)):
    f.write("group g%d %s\n"%(i*2+1,groups[i][0]))
    f.write("group g%d %s\n"%(i*2+2,groups[i][1]))
    f.write("\n")

# move and run spring restraint over the numebr of bins
for i in range(bins):
    for j in range(len(groups)):
        r=r0[j]+i*dr[j]
        xyz=np.zeros(3)
        xyz[direction]=r 
        if(j==1):
            xyz[2]=8.06
        f.write('fix s%d g%d spring couple g%d %s %6.2f %6.2f %6.2f %d\n'%(j,2*j+1,2*j+2,k[j],xyz[0],xyz[1],xyz[2],0)) #TODO: make work for multiple
    f.write('run 500000\n')
    for j in range(len(groups)):
        f.write('unfix s%d\n'%(j))
    f.write("\n")
f.write("write_data outfile1.data\n")
f.close()


