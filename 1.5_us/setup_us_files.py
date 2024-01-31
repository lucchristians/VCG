import numpy as np

data = np.loadtxt("summary_distancesxyz.dat")
c_u = open('confs_used.dat','w')
#TODO: set up to be more robust and order based on lowers to highest always
def closest_point(data,x0):
    diff = np.abs(data - x0)
    
    return np.argmin(diff)
x0=0.59
dx=0.1
total_points=30
for i in range(total_points):
    print(closest_point(data[:,1],x0+dx*i))
    c_u.write('%d\n'%(closest_point(data[:,1],x0+dx*i)))
c_u.close()
#curr_points=1

#print(int(data[0,0]))

#for i in range(len(data)):
    
    #if(x0-dx>np.abs(data[i,1])): #TODO: set up based on closest
    #    curr_points+=1
    #    x0-=dx
    #    print(int(data[i,0]))
#    if(curr_points==total_points):
#        break
