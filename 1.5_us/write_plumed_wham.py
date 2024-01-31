import numpy as np
import sys
import os

path=sys.argv[1]

if(path[-1]!='/'):
    path+='/'
print(path)
var=sys.argv[2]
if(var[-1]=='x' or var[-1]=='y' or var[-1]=='z'):
    cv=var[:-1]+'.'+var[-1]
variables=[]
restraints=[]
rep=0 #TODO: change to 0 if d_0 finished 
while(os.path.exists(path+'d_%d/us.dat'%(rep))):
    rep+=1

wham=open('plumed-wham.dat','w')

for i in range(rep):
    #if(i==0): #TODO: change to 0 if d_0 finished 
    #    continue
    with open(path+'d_%d/us.dat'%(i)) as f:
        for line in f:
            if(len(line)==0):
                continue
            elif(line[:3].lower()=='res'):
                if(var==line[:line.index(':')]):
                    restraints.append(line[line.index(':')+2:])
            elif(line.split()[1]=='COM' or line.split()[1]=='DISTANCE'): 
                present=False
                for v in variables:
                    if(v==line):
                        present=True
                        break
                if( not present):
                    variables.append(line)
for v in variables:
    wham.write(v)
for r in restraints:          
    wham.write(r)
wham.write('PRINT ARG=*.bias FILE=biases.dat STRIDE=1\n')
wham.write('PRINT ARG=%s FILE=res.dat STRIDE=1\n'%('dist'+cv[3:]))
wham.close()
