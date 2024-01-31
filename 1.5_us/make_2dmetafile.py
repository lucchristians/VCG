import numpy as np
import sys
import os

def read_coords(fileName,var): #only need distance along var
    data=[]
    with open(fileName) as f:
        #determine var fields to pull from
        first_line=True
        var_ind=[]
        res_ind=[]
        for line in f:
            dat=line.split()
            if(len(dat)==0):
                continue
            if(first_line):
                for i in range(len(var)):
                    res='res'+var[i][4:].replace('.','')+'.bias'
                    for j in range(len(dat)):
                        if(var[i]==dat[j]):
                            var_ind.append(j-2)
                        if(res==dat[j]):
                            res_ind.append(j-2)
                first_line=False
                continue
            #print(var_ind)
            data.append([float(dat[0]),float(dat[var_ind[0]]),float(dat[var_ind[1]])]) #TODO: find format
        data=np.array(data)
        return data

def find_window(fileName,var):
    restraints=[]
    info=[]
    for i in range(len(var)):
        res='res'+var[i][4:].replace('.','')
        with open(fileName) as f:
            for line in f:
                if(len(line)==0):
                    continue
                elif(line[:3].lower()=='res'):
                    if(res==line[:line.index(':')]):
                        dat=line.split()
                        distance=float(dat[3][dat[3].index('=')+1:])
                        spring_const=float(dat[4][dat[4].index('=')+1:])
                        info.append([spring_const,distance])
                        break
                        #restraints.append([spring_cosnt,distance])

                #elif(line.split()[1]=='COM' or line.split()[1]=='DISTANCE'):
                    #present=False
                    #for v in variables:
                    #    if(v==line):
                    #        present=True
                    #        break
                    #if( not present):
                    #    variables.append(line)
    return np.array(info)

rxn_coord1=sys.argv[1]
rxn_coord2=sys.argv[2]
rxn_coord=[rxn_coord1,rxn_coord2]
temperature=sys.argv[3]

#determine paths
curr_dir=os.getcwd()
if(curr_dir[-1]!='/'):
    curr_dir+='/'
num=0
while(os.path.exists(curr_dir+'d_%d'%(num))):
    num+=1

#
meta_file=open('metafile_2d.dat','w')
meta_file.write('# path r0 spring [correlation time] [temperature]\n')
for i in range(num):
    if(os.path.exists(curr_dir+'d_neg%d'%(i))):
        path=curr_dir+'d_neg%d/'%(i)
        timeseries=np.copy(read_coords(path+'colvar',rxn_coord))
        np.savetxt(path+'timeseries.dat',timeseries,fmt='%9.5f')
        window=np.copy(find_window(path+'us.dat',rxn_coord))
        meta_file.write("%s  %8.3f  %8.3f  %8.3f  %8.3f\n"%(path+'timeseries.dat',window[0,1],window[1,1],window[0,0],window[1,0]))
    path=curr_dir+'d_%d/'%(i)
    timeseries=np.copy(read_coords(path+'colvar',rxn_coord))
    np.savetxt(path+'timeseries.dat',timeseries,fmt='%9.5f')
    window=np.copy(find_window(path+'us.dat',rxn_coord))
    meta_file.write("%s  %8.3f  %8.3f  %8.3f  %8.3f\n"%(path+'timeseries.dat',window[0,1],window[1,1],window[0,0],window[1,0]))
    if(os.path.exists(curr_dir+'d_%d_5'%(i))):
        path=curr_dir+'d_%d_5/'%(i)
        timeseries=np.copy(read_coords(path+'colvar',rxn_coord))
        np.savetxt(path+'timeseries.dat',timeseries,fmt='%9.5f')
        window=np.copy(find_window(path+'us.dat',rxn_coord))
        meta_file.write("%s  %8.3f  %8.3f  %8.3f  %8.3f\n"%(path+'timeseries.dat',window[0,1],window[1,1],window[0,0],window[1,0]))
meta_file.close
