# -*- coding:utf-8 -*-
import os
import numpy as np
import sys
def alter(file, old_str, new_str):
	f = open(file,"r")
	file_data = ""
	for line in f:
	    if old_str in line:
	    	line = line.replace(old_str,new_str)
	    file_data += line
	f = open(file,"w")
	f.write(file_data)


prefix_name = sys.argv[1]
old_str = sys.argv[2]
filen = int(sys.argv[3])+1
new_str = ''

file_num = np.arange(0,filen)
for i in file_num:
	file = prefix_name+str(i).zfill(4)+'.xmf'
	print 'altering ' + file
	alter(file,old_str,new_str)