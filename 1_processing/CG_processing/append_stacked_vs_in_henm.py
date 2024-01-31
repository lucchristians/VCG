import numpy as np
import sys
import operator

base_result_str=sys.argv[1]
stack_result_str=sys.argv[2]
path=sys.argv[3]
if(path[-1]!='/'):
    path+='/'
#load in base and stacked result model
result=np.loadtxt(base_result_str,skiprows=1)
stack=np.loadtxt(stack_result_str)

#initialize results
r = open(path+"result_stacked.txt",'w')
r.write("    Atom I     Atom J         R0          K    Fluc_Ref. Fluc_Matched\n")

#append results into one array
results_new=[]
for resul in result:
    results_new.append(resul[:4])
for stac in stack:
    results_new.append(stac[[0,1,3,2]])
results_new=np.array(results_new)
#sort array based on ind i then index j
list1 =sorted(results_new, key=operator.itemgetter(0,1))

#print in a reasonable format
for i in range(len(list1)):
    results_=list1[i]
    r.write("%10d%11d%11.3f%11.3f%13.3f%13.3f\n"%(results_[0],results_[1],results_[2],results_[3],0.00,0.00))
r.close()
