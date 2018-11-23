import re
import numpy as np
import matplotlib.pyplot as plt
name = '/home/pzhang/chen/move-bed/'
name = name +'periodic.out'
f = open(name)
p = re.compile(r'\[.*?\]')
p1 = re.compile(r'-?\d+\.?\d*e?-\d\d')
num = []
for line in f.readlines():
    temp = re.findall(p,line)
    tt = []
    for t in temp:
        for r in re.findall(p1,t):
            tt.append(float(r))
    num.append(tt)

nnum = np.array(num)
for i in range(1,len(nnum)):
    ttemp = num[i]
    print ttemp[0]
    plt.plot(i,ttemp[0],'*')
plt.show()
