import numpy as np
import h5py as h5
import matplotlib.pyplot as plt 
import sys
from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("getZ.pyx"))
from getZ import *



prefix = "/media/pzhang/Elements/move-bed-tmp/macondo/"
prefix = prefix + "5e3Re_27.5Ga_0.3gap_a/"
prefix = prefix + "test_mvbed_c_"

Py = 21
Px = 160
R = 10
Np = Py*Px
print(Np)
num = np.arange(500,924 +1)

name = prefix + str(0).zfill(4) + ".h5"
f = h5.File(name)
Nx = np.array(f['Nx'])
Ny = np.array(f['Ny'])
Zz = np.zeros((Ny[0],Nx[0]))

for i in range(len(num)):
	name = prefix + str(num[i]).zfill(4) + ".h5"
	f = h5.File(name)
	# print("Nx ", Nx[0])
	# print("Ny ", Ny[0])
	x = np.arange(Nx[0])
	y = np.arange(Ny[0])
	xx , yy= np.meshgrid(x,y,indexing='xy')
	pos = np.array(f['Pposition'])

	px = np.array(pos[0:-2:3])
	py = np.array(pos[1:-1:3])
	# for ip in range(Np):
		# Z[np.where(((xx-px[ip])**2+(yy-py[ip])**2)<=R*R)]=1
	Z = getZ(Nx,Ny,Np,R,px,py)
	Zz = Zz+Z
	print(i)
Zz = Zz/len(num)
np.save("Zz.npy",Zz)
# fig,ax = plt.subplots(1,1)
# c = ax.pcolor(Zz)
# fig.colorbar(c,ax=ax)
# plt.show()



