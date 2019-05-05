import os
import sys

dir_name = sys.argv[1]
file_prefix = sys.argv[2]

first = int(sys.argv[3])
second = int(sys.argv[4])

for i in range(first,second+1):
	nameh5 = dir_name + file_prefix + str(i).zfill(4) + ".h5"
	print "prepare to rm " + nameh5
	os.remove(nameh5)
	namexmf = dir_name + file_prefix + str(i).zfill(4) + ".xmf"
	print "prepare to rm " + namexmf
	os.remove(namexmf)

