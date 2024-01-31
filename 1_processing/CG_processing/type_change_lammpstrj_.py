import numpy as np
from mscg import *
import sys

input_file=sys.argv[1]
output_file=sys.argv[2]

traj=Trajectory(input_file, fmt='lammpstrj')
trajout=Trajectory(output_file, fmt='lammpstrj', mode='w')
while(traj.read_frame()):
    trajout.natoms=traj.natoms + 5
    trajout.box=traj.box
    #define type list
    trajout.t=np.arange(0,len(traj.t))%54+1
    #add in stacked vs sites (on bottom chains)
    trajout.x=traj.x
    trajout.write_frame()
    #pos=traj.x
    #pos_array.append(np.copy(pos))i


