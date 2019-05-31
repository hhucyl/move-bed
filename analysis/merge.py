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
def process(fn):
	ffn = ((fn.split('/'))[-1].split('_'))[-1].split('.')
	fffn = ffn[0]
	return int(fffn)
	


	
prefix_name = sys.argv[1]
print "The directory is: " + prefix_name
n = len(sys.argv)
file_name = sys.argv[2]
print "Merge all file to this name: " + file_name
file_list = os.listdir(prefix_name)
# print file_list
nf = np.zeros(n-2,dtype=int)
flist = []
for ii in np.arange(2,n):
	ftemp = []	
	for i in range(len(file_list)):
		if(re.search(sys.argv[ii] + '\d+',file_list[i])):
			nf[ii-2] = nf[ii-2]+1
			fftemp = [process(file_list[i]), file_list[i], len(file_list[i])+process(file_list[i])]
			ftemp.append(fftemp)
	ftemp = sorted(ftemp,key=lambda x: x[0])#identify xmf h5 maybe not in order!!!!
	ftemp = sorted(ftemp,key=lambda x: x[2])

	flist.append(ftemp)
nf = nf/2
print nf
# print flist[1]
print flist[0][-1][0]
print len(flist[0])/2-1
nnn = len(flist[0])/2-1
if flist[0][0][0]>0:

	for i in np.arange(nnn+1):
	
		h5name = prefix_name + flist[0][2*i][1]	
		h5name1 = prefix_name + file_name + str(i).zfill(4) + ".h5"
		print "Renaming " + h5name + " to " + h5name1
		flist[0][2*i] = [process(h5name1),h5name1.split('/')[-1]]
		xmfname = prefix_name + flist[0][2*i+1][1]		
		xmfname1 = prefix_name + file_name + str(i).zfill(4) + ".xmf"
		print "Renaming " + xmfname + " to " + xmfname1
		flist[0][2*i+1] = [process(xmfname1),xmfname1.split('/')[-1]]
		# os.rename(h5name,h5name1)		
		# os.rename(xmfname,xmfname1)
		# alter(xmfname1,xmfname)
nnn = flist[0][-1][0]	
for i in np.arange(3,n):
	for ii in range(len(flist[i-2])/2):
		nnn = nnn+1
		#h5name = prefix_name + sys.argv[i] + str(ii).zfill(4) + ".h5"
		h5name = prefix_name + flist[i-2][2*ii][1]		
		h5name1 = prefix_name + file_name + str(nnn).zfill(4) + ".h5"	
		print "Renaming " + h5name + " to " + h5name1
		
		# xmfname = prefix_name + sys.argv[i] + str(ii).zfill(4) + ".xmf"
		xmfname = prefix_name + flist[i-2][2*ii+1][1]		
		xmfname1 = prefix_name + file_name + str(nnn).zfill(4) + ".xmf"
		print "Renaming " + xmfname + " to " + xmfname1	
		os.rename(h5name,h5name1)		
		os.rename(xmfname,xmfname1)
		alter(xmfname1,xmfname)



