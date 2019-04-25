import os
import sys
import numpy as np
import re

def alter(xf1,xf):
	f = open(xf1,"r")
	file_data = ""
	xff = (xf.split('/')[-1]).split('.')[0]
	xff1 = (xf1.split('/')[-1]).split('.')[0]
	
	for line in f:
		if xff in line:
			# pos = line.index(file_name)
			# line = line[pos:]
			# print xff
			
			line = line.replace(xff,xff1)
		file_data += line
	f = open(xf1,"w")
	f.write(file_data)
	f.close()

prefix_name = sys.argv[1]
print "The directory is: " + prefix_name
n = len(sys.argv)
file_name = sys.argv[2]
print "Merge all file to this name: " + file_name
file_list = os.listdir(prefix_name)
# print file_list
nf = np.zeros(n-2,dtype=int)
for ii in np.arange(2,n):
	for i in range(len(file_list)):
		if(re.search(sys.argv[ii] + '\d+',file_list[i])):
			nf[ii-2] = nf[ii-2]+1
nf = nf/2
print nf
nn = nf[0]
for i in np.arange(3,n):
	for ii in range(nf[i-2]):
		h5name = prefix_name + sys.argv[i] + str(ii).zfill(4) + ".h5"
		h5name1 = prefix_name + file_name + str(nn).zfill(4) + ".h5"
		print "Renaming " + h5name + " to " + h5name1
		os.rename(h5name,h5name1)
		xmfname = prefix_name + sys.argv[i] + str(ii).zfill(4) + ".xmf"
		xmfname1 = prefix_name + file_name + str(nn).zfill(4) + ".xmf"
		print "Renaming " + xmfname + " to " + xmfname1	
		os.rename(xmfname,xmfname1)

		alter(xmfname1,xmfname)
		nn = nn + 1

