import numpy as np
import sys

base_result_str=sys.argv[1]
stack_result_str=sys.argv[2]
path=sys.argv[3]
#load in base and stacked result model
result=np.loadtxt(base_result_str,skiprows=1)
stack=np.loadtxt(stack_result_str)

#initialize results
r = open(path+"result_r0.txt",'w')
r.write("    Atom I     Atom J         R0          K    Fluc_Ref. Fluc_Matched\n")

#append results into one array
results_new=[]
for resul in result:
    in_const=False
    for i in range(len(stack)):
        if(np.all(resul[:2]==stack[i,:2])):
            print(i)
            in_const=True
            ind=i
            break
    if(in_const):
        temp=stack[ind,:4]
        temp_=np.copy(stack[ind,3])
        temp[3]=np.copy(stack[ind,2])
        temp[2]=temp_
        results_new.append(temp)
        
    else:
        results_new.append(resul[:4])
results_new=np.array(results_new)
#sort array based on ind i then index j
#list1 =sorted(results_new, key=operator.itemgetter(1,2))

#print in a reasonable format
for i in range(len(results_new)):
    results_=results_new[i]
    r.write("%10d%11d%11.3f%11.3f%13.3f%13.3f\n"%(results_[0],results_[1],results_[2],results_[3],result[i,4],result[i,5]))
r.close()
