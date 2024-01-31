import numpy as np
import sys

outputName=sys.argv[1]
windows=int(sys.argv[2])
rinit=float(sys.argv[3])
rfin=float(sys.argv[4])
row=sys.argv[5]
dt = 30 #fs
trans_steps = 50000 #n steps to move
equil_steps = 450000 #n steps to stay in place, discard some of this data too

dr = (rfin - rinit) / (windows-1)

f = open(outputName,'w')

# define group that doesn't move
f.write("fix vsi res move linear 0 0 0\n")

# define moving groups
for i in range(windows):
    vel = dr / (trans_steps * dt)
    if(row=='x'):
        f.write(f"fix d1 vsj move linear {vel} 0.0 0.0\n")
    elif(row=='y'):
        f.write(f"fix d1 vsj move linear 0.0 {vel} 0.0\n") 
    elif(row=='z'):
        f.write(f"fix d1 vsj move linear 0.0 0.0 {vel}\n")
    f.write(f'run {trans_steps}\n')
    f.write('unfix d1\n')

    f.write(f"fix d1 vsj move linear 0.0 0.0 0.0\n")
    f.write(f'run {equil_steps}\n')

f.write("write_data outfile1.data\n")
f.close()


