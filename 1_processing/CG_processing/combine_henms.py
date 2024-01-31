import numpy as np
res=[]

flex_region=np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,54])+1
imflex_region=np.arange(17,59,dtype=int)+1
#vs_f_f=np.array([22])
#vs_t_f=np.array([55])
#vs_f_if=np.array([38,39,40,41,42])
#vs_t_if=np.array([55,56,57,58,59])

def replace_vs_sites(data_,vs_f):
    data=np.copy(data_)
    for i in reversed(range(len(vs_f))):
        for j in range(2):
            mask = data[:,j]==i+1
            data[mask,j]=vs_f[i]
    return data

new_henm=open("rcut_combine/result.txt",'w')
new_henm.write("    Atom I     Atom J         R0          K    Fluc_Ref. Fluc_Matched\n")
data_if = np.loadtxt("rcut_15imflex/result.txt",skiprows=1)
data_if = np.copy(data_if)
data_f = np.loadtxt("rcut_15flex/result.txt",skiprows=1)
data_f = np.copy(data_f)

#length_f=len(data_f)-len(vs_t_f)
#map vs
#for flex
df = replace_vs_sites(data_f,flex_region)

# shift indexing
#data_if[:,0]=data_if[:,0]+
#data_if[:,:2]=data_if[:,:2]+length_f
#vs_f_if+=length_f
#for imflex
dif = replace_vs_sites(data_if,imflex_region)

#print(dif[:,:2])
#print(df[:,:2])
# add all flex first
start_if=0
indexes_used=[]
data=[]
temp_f=[]
for k in range(len(df)):
    averaged=False
    for l in range(len(dif)):
        #print(df[k,0],dif[l,0],df[k,1],dif[l,1])
        
        if(df[k,0]==dif[l,0] and df[k,1]==dif[l,1]):
            frame=df[k,:]
            frame[2:]=(df[k,2:]+dif[l,2:])/2
            temp_f.append(np.copy(frame))
            averaged=True
            indexes_used.append(l)
            print(df[k,:2])
            start_if+=1
            break 
    if(not averaged):
        #print(df[k,:2])
        data.append(df[k,:])
#copy in imflexable trajectories
for k in range(len(dif)):
    used=False
    for l in range(len(indexes_used)):
        if(indexes_used[l]==k):
            used=True
            data.append(temp_f[l])
            break
    if(used):
        continue
    data.append(dif[k,:])    
#data.append(data_f)
#print(indexes_used)
#print(np.shape(data))
data=np.array(data)
#exit()
for k in range(len(data)):
    #if((data[k,5]-data[k,4])**2<0.1):
    new_henm.write("%10d %10d %10.3f %10.3f   %10.3f   %10.3f\n"%(data[k,0],data[k,1],data[k,2],data[k,3],data[k,4],data[k,5]))
new_henm.close()
#print(res)
