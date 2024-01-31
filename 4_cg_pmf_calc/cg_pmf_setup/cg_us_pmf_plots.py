import numpy as np
import sys
import matplotlib.pyplot as plt

path=sys.argv[1]
if(path[-1]!='/'):
    path+='/'
if(path=='vs/inter/'):
    print("modified")
    special=True #had to rerun which changed log format
else:
    special=False
logFile=sys.argv[2]
outputName=sys.argv[3]
r0=float(sys.argv[4])
rf=float(sys.argv[5])
def find_potential(path,logFile):
    distances=[]
    potentials=[]
    skip_first_step=True # is minimize step
    gather_data=False
    print(path+logFile)
    with open(path+logFile) as f:
        for line in f:
            data=line.split()
            if(len(data)==0):
                continue
            if(data[0]=='Step'):
                if(skip_first_step):
                    skip_first_step=False
                    continue
                gather_data=True
                #initalize frame
                potential=[]
                continue
            if(gather_data):
                if(data[0]=='Loop'): #process data
                    #print(np.shape(potential))
                    potential=np.array(potential)
                    if(not special):
                        potentials.append(np.copy(potential[250:])) # ignore movement time
                    else:    
                        potentials.append(np.copy(potential)) # ignore movement time
                    
                    gather_data=False
                    continue 
                #print(path+logFile,data[4])
                potential.append(float(data[4]))
    if(special):
        potentials=np.array(potentials).reshape(30,-1)[:,250:]
    else:
        potentials=np.array(potentials)
    
    print(logFile,np.shape(potentials))
    np.savetxt(path+logFile+'_potential.dat',potentials)
    return np.copy(potentials)

pots=np.zeros([30,0])
for i in range(4):
    if(special and i>2):
        continue
    temp = find_potential(path,logFile+'_%d'%(i+1))
    pots=np.hstack([pots,temp])
    print(np.shape(pots))
    print("")
r=np.linspace(r0,rf,num=30)
pot_mean=np.mean(pots,axis=1)

#shift potential to have end at zero
pot_mean=pot_mean-pot_mean[-1]
pot_std=np.std(pots,axis=1)
np.savetxt(outputName[:-4]+'.dat',np.array([r,pot_mean,pot_std]))
plt.figure()
plt.plot(r,pot_mean,'-k')
plt.plot(r,pot_mean-pot_std,'--k')
plt.plot(r,pot_mean+pot_std,'--k')
#plt.show()
plt.savefig(outputName)


