# -*- coding:utf-8 -*-
import os
import numpy as np
import sys
# def alter(file, old_str, new_str):
# 	f = open(file,"r")
# 	file_data = ""
# 	for line in f:
# 	    if old_str in line:
# 	    	line = line.replace(old_str,new_str)
# 	    file_data += line
# 	f = open(file,"w")
# 	f.write(file_data)

def alter(file,file_name):
	f = open(file,"r")
	file_data = ""
	for line in f:
		if file_name in line:
			pos = line.index(file_name)
			line = line[pos:]
		file_data += line
	f = open(file,"w")
	f.write(file_data)

prefix_name = sys.argv[1]
# old_str = sys.argv[2]
fistart = int(sys.argv[2])
fiend = int(sys.argv[3])+1
# new_str = ''
file_name = prefix_name.split('/')[-1]
# print file_name
# print fistart
# print fiend
file_num = np.arange(fistart,fiend)
for i in file_num:
	file = prefix_name+str(i).zfill(4)+'.xmf'
	print 'altering ' + file
	# alter(file,old_str,new_str)
	alter(file,file_name)