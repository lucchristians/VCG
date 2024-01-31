import numpy as np 
import sys
#defines the 25, 27, and 29th residues to be restrained
def restraint_residues(groFilePath):
    goal=[25,27,29]
    res_inds=[]
    with open(groFilePath) as f:
        count=0
        for line in f:
            if(count<2):
                count+=1
                continue
            
            resid=int(line[0:5])
            if(resid>50):
                break
            atom=str(line[13:15])
            index=int(line[15:20])
            
            for ind in goal:
                
                if(resid==ind and atom=='CA'):
                    res_inds.append(index)
                    break
            count+=1
    return np.array(res_inds)
            
#inputs: python script.py path_to_topol_files N_topol_files char1 char2 ..  charN
path=sys.argv[1]
if(path[-1]!='/'):
    path+='/'
topols_input=[]
topols_output=[]
# read in characters
for i in range(int(sys.argv[2])):
    # topols_input.append(path+"/posre_protein_chain_%s.itp")
    topols_output.append(path+"posre_Protein_chain_%s_lock.itp"%(sys.argv[3+i]))
res_ind = restraint_residues(path+"minimize.gro")
print(res_ind)
for i, input_file in enumerate(topols_output):
    print(topols_output[i])
    o = open(topols_output[i],'w')
    # write the header comments and section definition
    o.write("; In this topology include file, you will find position restraint \n")
    o.write("; entries for all the heavy atoms in your original pdb file. \n")
    o.write("; This means that all the protons which were added by pdb2gmx are \n")
    o.write("; not restrained. \n")
    o.write(" \n")
    o.write("[ position_restraints ]\n")
    o.write("; atom  type      fx      fy      fz\n")
    
    # write the atoms that need to be restrained
    for res in res_ind:
        o.write("%6d%6d%6d%6d%6d\n"%(res,1,1000,1000,1000))
    o.close()
    
    # append include statement to end of topol for the specific chain
    #check if the specific include statement is at teh end of the file already
    already_included=False
    char=sys.argv[3+i]
    include_format=["#ifdef POSRESL\n",'#include "posre_Protein_chain_%s_lock.itp"\n'%(char),"#endif\n"]
    with open(path+"topol_Protein_chain_%s.itp"%(char)) as f:
        lines_in_a_row=0
        max_lines=0
        for line in f:
            present=False
            for string in include_format:
                if(string==line):
                    present=True
                    break
            if(present):
                lines_in_a_row+=1
            else:
                if(lines_in_a_row>max_lines):
                    max_lines=lines_in_a_row
                lines_in_a_row=0
        if(lines_in_a_row==3):
            already_included=True

    #if not add it/ if so ignore it
    if(not already_included):
        f = open(path+"topol_Protein_chain_%s.itp"%(char),'a')
        f.write(include_format[0])
        f.write(include_format[1])
        f.write(include_format[2])
        f.close()
