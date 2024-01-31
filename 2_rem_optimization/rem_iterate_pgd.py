#!/usr/bin/env python

import numpy as np
from mscg import *
import matplotlib.pyplot as plt
import os 

# REDO of REM method
# Explicitly calculate expectation of AA/CG using mean/var
# Will need to reprocess the AA/CG trajectory each time

### 
#   Global Variables
###

# Cutoff for calculating for histogram in Angstroms
binwidth=0.05 #A
max=20 #20.0
bins=int(max/binwidth)

def create_hist(R,data,n_frames,Np):
    shape=np.shape(data)[0]
    rlo=0.0
    rhi=max
    g , x = np.histogram(data,bins=bins,range=(rlo,rhi))
    return (x[:-1]+x[1:])/2.0, g/n_frames/Np

# Physical constants
Temp=300 #K, this was wrong in Ethan's script
kb=1.987204E-3 #kcal/mol/K
beta=1.0/(kb*Temp) #mol/kcal

# Max change of a varible in one iteration
max_change_A=0.50  
max_change_B=0.02  
min_A = 0.01
min_B = 0.1 #if B gets smaller than this, you are exploring an unphysically long-range interaction

noise_width_A = 0.050 #random noise added to each update, I set this to 10% of the max update
noise_width_B = 0.002

###
#   Functions
###

# Defining derivative
dudA = lambda K,A,B,x,R: (-np.exp(-B*(x-R)**2))
dudB = lambda K,A,B,x,R: A*((x-R)**2)*np.exp(-B*(x-R)**2) 

du2dA2 = lambda K,A,B,x,R: np.zeros_like(x)
du2dB2 = lambda K,A,B,x,R: -A*((x-R)**4)*np.exp(-B*(x-R)**2)


###
#   REM Script
###

# Check for most recent REM iteration (regardless of if it is finished or not)
num=1
path = os.getcwd() 
K_const_ind=np.loadtxt("K_const_ind.txt")
while(os.path.exists(path+"/rem_iteration"+str(num)+"/outfile.in.settings")):
    num+=1
num-=1

# learning rate schedule
if (num<=200):
    # change A and B
    chi_A = 0.1
    chi_B = 0.1
    chi_K = 0.1
else :
    # change both to lower rate to help find local minima (your mileage may vary)
    chi_A = 0.01
    chi_B = 0.01
    chi_K = 0.01

# Gets initial constants from last iteration
lambda_k1=[]
lambda_k2=[]
if (num>0):
    with open("rem_iteration"+str(num)+"/outfile.in.settings",'r') as old_settings:
        for line in old_settings:
            info = line.split()
            
            if(len(info)<4): # Filters empty lines
                continue
            
            if(info[3] == "gauss"):
                lambda_k1.append([float(info[4]),float(info[5])])
            
            if(info[0] == 'bond_coeff'):
                # Finds location of the vs-associated bond
                is_a_K_const=False 
                for k in range(len(K_const_ind)): 
                    if(int(info[1])==int(K_const_ind[k])):
                        is_a_K_const=True
                        break
                
                if(is_a_K_const):
                    lambda_k2.append(info[2])
else:
    print("script can not be used without an iteration")
    exit()

lambda_k=np.hstack([np.reshape(lambda_k2,(len(lambda_k2),1)),np.array(lambda_k1)]).astype(float)


# Load in, reads, & process trajectories
cg_traj = Trajectory("rem_iteration"+str(num)+"/sub1.lammpstrj", fmt='lammpstrj')
cg_pos_array = []
cg_cell_lens = cg_traj.box
while(cg_traj.read_frame()):
    atomtype_cg=cg_traj.t
    pos=cg_traj.x
    cg_pos_array.append(np.copy(pos))
cg_n_frames=len(cg_pos_array)
cg_pos_array=np.array(cg_pos_array)

aa_traj = Trajectory("reference_AA/2hexamers_combine.lammpstrj", fmt='lammpstrj')
aa_pos_array = []
aa_cell_lens = aa_traj.box
while(aa_traj.read_frame()):
    atomtype_aa=aa_traj.t
    pos=aa_traj.x
    aa_pos_array.append(np.copy(pos))
aa_n_frames=len(aa_pos_array)
aa_pos_array=np.array(aa_pos_array)

# Indexing for gauss fitting
txt=np.loadtxt('pairhist.txt')
if(len(np.shape(txt))==1):
    txt = txt.reshape(1,2)

# Shift for distributions
const=np.loadtxt('reference_AA/const_AA.txt')

# Check if enough pairs are specified to parameterize all gaussians
if(np.shape(txt)[0] != np.shape(lambda_k)[0]):
    Exception("you don't have enough pairs to parameterize all constants: (pairs: %d | parameters: %d"%(np.shape(txt)[0],np.shape(lambda_k)[0]))
    exit()

# Create/append output information 
if(num == 1): 
    loss = open("logs/loss.txt",'w')
    loss_mse = open("logs/loss_mse.txt",'w')
    constants = open("logs/constants.txt",'w')
    debug_log_A = open("logs/debug_log_A.txt", "w")
    debug_log_B = open("logs/debug_log_B.txt", "w")
    
    # Create constants file header
    constants.write("#iter ")
    for i in range(len(txt)):
        constants.write("K%d A%d B%d "%(i,i,i))
    constants.write("\n")
    loss.write("#iter Srel\n") 
else:
    loss = open("logs/loss.txt",'a')
    loss_mse = open("logs/loss_mse.txt",'a')
    constants = open("logs/constants.txt",'a')
    debug_log_A = open("logs/debug_log_A.txt", "a")
    debug_log_B = open("logs/debug_log_B.txt", "a")

debug_log_B.write("--------------- Iteration: %d --------------- \n" % num)
debug_log_A.write("--------------- Iteration: %d --------------- \n" % num)


# REM update
probs_AA = np.loadtxt("reference_AA/AA_hist.txt")
ranges=[]
probs=[]
new_const = []
for i in range(len(txt)):
    # Shape array from [frames, number of types, 3] to [frames * number of types, 3]
    atom1=txt[i,0]
    atom2=txt[i,1]

    #current forcefield constants
    K_curr = lambda_k[i][0]
    A_curr = lambda_k[i][1]
    B_curr = lambda_k[i][2]

    aa_duda = 0.0
    aa_du2da2 = 0.0
    aa_dudb = 0.0
    aa_du2db2 = 0.0

    b = aa_pos_array[:,atomtype_aa==atom1,:]
    c = aa_pos_array[:,atomtype_aa==atom2,:]

    for j in range(np.shape(b)[1]):
        for k in range(np.shape(c)[1]):

            if(j==k):
                continue
            bc = b[:,j,:] - c[:,k,:]
            dr = bc - aa_cell_lens * np.rint( bc / aa_cell_lens )
            rr = np.linalg.norm(dr, axis=1) # (nframes,1)

            aa_duda += dudA( K_curr, A_curr, B_curr, rr, const[i] )
            aa_du2da2 += du2dA2( K_curr, A_curr, B_curr, rr, const[i] )

            aa_dudb += dudB( K_curr, A_curr, B_curr, rr, const[i] )
            aa_du2db2 += du2dB2( K_curr, A_curr, B_curr, rr, const[i] )
            

    cg_duda = 0.0
    cg_du2da2 = 0.0
    cg_dudb = 0.0
    cg_du2db2 = 0.0

    Np = 0
    b = cg_pos_array[:,atomtype_cg==atom1,:]
    c = cg_pos_array[:,atomtype_cg==atom2,:]
    d = []

    for j in range(np.shape(b)[1]):
        for k in range(np.shape(c)[1]):

            if(j==k):
                continue
            bc = b[:,j,:]-c[:,k,:]
            dr = bc - cg_cell_lens * np.rint( bc / cg_cell_lens )
            d.append(dr)
            rr = np.linalg.norm(dr, axis=1) # (nframes,1)

            cg_duda += dudA( K_curr, A_curr, B_curr, rr, const[i] )
            cg_du2da2 += du2dA2( K_curr, A_curr, B_curr, rr, const[i] )

            cg_dudb += dudB( K_curr, A_curr, B_curr, rr, const[i] )
            cg_du2db2 += du2dB2( K_curr, A_curr, B_curr, rr, const[i] )
            
            Np+=1

    d = np.swapaxes(np.array(d),0,1)
    e = np.linalg.norm(d,axis=2)
    f = e.flatten()
    x, h = create_hist(const[i],f,cg_n_frames,Np)

    # because of the way counting is done, these are overcounted by 2x
    if(atom1 == atom2) :
        cg_duda /= 2.0
        cg_du2da2 /= 2.0
        cg_dudb /= 2.0
        cg_du2db2 /= 2.0

        aa_duda /= 2.0
        aa_du2da2 /= 2.0
        aa_dudb /= 2.0
        aa_du2db2 /= 2.0

    # now perform rem update

    # Determines the gradient and hessian for the gaussian magnitude
    dUdA_AA = aa_duda.mean()
    dUdA_CG = cg_duda.mean()
    d2UdA2_AA = aa_du2da2.mean()
    d2UdA2_CG = cg_du2da2.mean()
    dUdA_var_CG = cg_duda.var()
    
    grad_A = beta * (dUdA_AA - dUdA_CG)
    hessian_A = beta * (d2UdA2_AA - d2UdA2_CG + beta * dUdA_var_CG )
    debug_log_A.write("vs_type: %d\n"%(i))
    debug_log_A.write("Gradient: %f\n" % grad_A)
    debug_log_A.write("Gradient terms: %f %f\n" % (dUdA_AA, dUdA_CG) )
    debug_log_A.write("Hessian: %f\n" % hessian_A)
    debug_log_A.write("Hessian terms: %f %f %f\n" % (d2UdA2_AA, d2UdA2_CG, dUdA_var_CG))

    # only use PGD update
    debug_log_A.write("in perturbed gradient descent\n")
    lambda_A=lambda_k[i][1] - chi_A * grad_A + np.random.normal( loc=0.0, scale=noise_width_A )

    if (lambda_A < min_A ) :
        lambda_A = min_A
        
    # Calculates expectations for gaussian width derivatives
    # Determines the gradient and hessian for the gaussian width
    dUdB_AA = aa_dudb.mean()
    dUdB_CG = cg_dudb.mean()
    d2UdB2_AA = aa_du2db2.mean()
    d2UdB2_CG = cg_du2db2.mean()
    dUdB_var_CG = cg_dudb.var()

    grad_B = beta * (dUdB_AA - dUdB_CG)
    hessian_B = beta * (d2UdB2_AA - d2UdB2_CG + beta * dUdB_var_CG)
    debug_log_B.write("vs_type: %d\n"%(i))
    debug_log_B.write("Gradient: %f\n" % grad_B)
    debug_log_B.write("Gradient terms: %f %f\n" % (dUdB_AA, dUdB_CG) )
    debug_log_B.write("Hessian: %f\n" % hessian_B)
    debug_log_B.write("Hessian terms: %f %f %f\n" % (d2UdB2_AA, d2UdB2_CG, dUdB_var_CG))

    debug_log_B.write("in perturbed gradient descent\n")
    lambda_B=lambda_k[i][2] - chi_B * grad_B + np.random.normal( loc=0.0, scale=noise_width_B )
    if (lambda_B < min_B ) :
        lambda_B = min_B

    if(np.abs(lambda_A-lambda_k[i][1])>max_change_A):
        lambda_A=lambda_k[i][1]+(lambda_A-lambda_k[i][1])/(np.abs(lambda_A-lambda_k[i][1]))*max_change_A
    if(np.abs(lambda_B-lambda_k[i][2])>max_change_B):
        lambda_B=lambda_k[i][2]+(lambda_B-lambda_k[i][2])/(np.abs(lambda_B-lambda_k[i][2]))*max_change_B
    lambda_K=lambda_A*lambda_B
    new_const.append([lambda_K,lambda_A,lambda_B])
    
    ranges.append(np.copy(x))
    probs.append(np.copy(h))

ranges=np.array(ranges)
probs=np.array(probs)

np.savetxt(f"rem_iteration{num}/CG_hist.txt", probs)
np.savetxt(f"rem_iteration{num}/CG_r.txt", ranges)
    
# Update simulation constants
constants.write("%d "%(num))
for i in range(np.shape(new_const)[0]):
    for j in range(np.shape(new_const)[1]):
        print("%7.4f"%(new_const[i][j]))
        constants.write("%7.4f "%(new_const[i][j]))
constants.write("\n")
constants.close()

# Calculates loss
loss.write("%d "%(num))
loss.write("%6.4f "%(chi_A))
loss.write("%6.4f "%(chi_B))
tot=0
for i in range(len(probs)):
    mask= np.all(np.concatenate(((probs[i]>0.0).reshape(-1,1),(probs_AA[i]>0.0).reshape(-1,1)),axis=1),axis=1)
    Srel = np.sum(probs_AA[i,mask]*np.log(probs_AA[i,mask]/probs[i,mask])) 
    loss.write("%10.5f "%(Srel))
    tot+=Srel
loss.write("%10.5f"%(tot))
loss.write("\n")
loss.close()

# Calculate mse loss
loss_mse.write("%d "%(num))
loss_mse.write("%6.4f "%(chi_A))
loss_mse.write("%6.4f "%(chi_B))
tot=0
for i in range(len(probs)):
    diff_p = probs_AA[i] - probs[i]
    Srel = np.sum( diff_p**2 )
    loss_mse.write("%10.5f "%(Srel))
    tot+=Srel
loss_mse.write("%10.5f"%(tot))
loss_mse.write("\n")
loss_mse.close()

