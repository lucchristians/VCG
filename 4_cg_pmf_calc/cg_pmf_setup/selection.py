import numpy as np

dr=62
i=346
f=363
print('set sel%d [atomselect top "index '%(4),end='')
for j in range(5):
    print('%d to %d '%(i+dr*j,f+dr*j),end='')
print('"]\n')
print('set com%d [measure center $sel%d weight mass]'%(4,4))
