import numpy as np
cimport numpy as np
def getZ(int Nx, int Ny, int Np, int R, np.ndarray px, np.ndarray py):
	cdef np.ndarray Z = np.zeros([Ny,Nx],dtype=np.int)
	cdef int ip, ixs, ixe, iys, iye
	for ip in range(Np):
		ixs = max(np.floor(px[ip]-R-3),0)
		ixe = min(np.ceil(px[ip]+R+3),Nx)
		iys = max(np.floor(py[ip]-R-3),0)
		iye = min(np.ceil(py[ip]+R+3),Ny)
		for xx in np.arange(ixs,ixe):
			for yy in np.arange(iys,iye):
				if(((xx-px[ip])**2+(yy-py[ip])**2)<=R*R):
					Z[yy,xx] = 1
	# 	print(ip)
	# ip = 44
	# ixs = max(np.floor(px[ip]-R-3),0)
	# ixe = min(np.ceil(px[ip]+R+3),Nx)
	# iys = max(np.floor(py[ip]-R-3),0)
	# iye = min(np.ceil(py[ip]+R+3),Ny)
	# print(ixs)
	# print(ixe)
	# print(iys)
	# print(iye)
	# for xx in np.arange(ixs,ixe):
	# 	for yy in np.arange(iys,iye):
	# 		if(((xx-px[ip])**2+(yy-py[ip])**2)<=R*R):
	# 			print(xx," ",yy)
	# 			Z[xx,yy] = 1
	return Z

