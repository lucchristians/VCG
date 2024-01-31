import numpy as np
import sys
import os

windows=int(sys.argv[1])
rinit=float(sys.argv[2])
rfin=float(sys.argv[3])
rorig = float(sys.argv[4])
coord = sys.argv[5]
#potential files to copy
fileNames=['input.pull','outfile.data','system.init','outfile.in.settings','exclusion.ff','coul.ff','hydrophobic.ff','intra_groups']

k_const = 0.5 #kcal/mol/A^2 for 1/2*K*(r-r0)^2

dr = (rfin - rinit) / (windows-1)

# define windows
for i in range(windows):

    xx = rinit + dr*i
    disp = xx - rorig
    
    os.mkdir(f"window_{i}")
    f = open(f"window_{i}/spring.settings",'w')
    if(coord=='x'):
        f.write(f"displace_atoms m1 move  {disp}  0.0  0.0\n")
        f.write(f"fix s1 g1 spring couple g2  {k_const}  {xx}   NULL  NULL   0.0\n")
    elif(coord=='y'):
        f.write(f"displace_atoms m1 move  0.0  {disp}  0.0\n")
        f.write(f"fix s1 g1 spring couple g2  {k_const}  NULL  {xx}   NULL   0.0\n")
    elif(coord=='z'):
        f.write(f"displace_atoms m1 move  0.0  0.0  {disp}\n")
        f.write(f"fix s1 g1 spring couple g2  {k_const}  NULL  NULL   {xx}   0.0\n")
    else:
        print('coord not properly set; must be "x", "y", or "z"') 
    for j in range(len(fileNames)):
        name=fileNames[j]
        if(os.path.exists(f"{name}")):
            os.symlink(f"../{name}", f"window_{i}/{name}")
       
    
f.close()


