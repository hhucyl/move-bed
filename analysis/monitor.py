import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

prefix = "/PZ_Q/move-bed-tmp/"
prefix = prefix + "5e3Re_26.0Ga_0.3gap_a/"
prefix = prefix + "test_mvbed_c_1_"
Py = 21
Px = 160
Np = Py*Px
num = np.arange(89,922+1)
nn = []
Vx = []
Vbx = []
Vpx = []
for i in range(len(num)):
	name = prefix + str(num[i]).zfill(4) + ".h5"
	f = h5.File(name)
	Nx = np.array(f['Nx']);
	Ny = np.array(f['Ny']);
	vel = np.array(f['Velocity_0']);
	vx = vel[0:-2:3]
	vx = vx.reshape((Ny[0],Nx[0]))
	vvx = np.abs(np.average(vx,axis=1))
	
	plt.subplot(3,2,1)
	plt.plot(vvx,np.arange(Ny[0]))
	plt.xlabel(r'$u_x$')
	plt.ylabel(r'depth')
	# plt.hold(True)

	plt.subplot(3,2,2)
	nn.append(num[i])
	vx = np.array(vx)
	idx = int(Nx[0]/2)
	vvvx = vx[:,idx]
	Vx.append(np.abs(vvvx.sum()))
	plt.plot(nn,Vx,'*')
	plt.ylabel(r'$\sum u_x$')
	plt.xlabel(r'num')
	# plt.hold(False)

	
	ax = plt.subplot(3,2,3)
	pos = np.array(f['Pposition'])
	ppy = pos[1:3*Py*Px-1:3].reshape(Py,Px)
	py = np.max(pos[1:-1:3])
	plt.plot(vvx[:int(py)],np.arange(np.size(vvx[:int(py)])))
	# ax.set_ylim(np.min(vvx[:int(py)]),np.max(vvx[:int(py)]))
	plt.xlabel(r'$u_x^{bed}$')
	plt.ylabel(r'depth')
	# plt.hold(True)
	
	ax = plt.subplot(3,2,4)
	Vbx.append(vvx[:int(py)].sum())
	plt.plot(nn,Vbx,'*')
	plt.ylabel(r'$\sum u_x^{bed}$')
	plt.xlabel(r'num')
	# plt.hold(False)
	
	plt.subplot(3,2,5)
	pv = np.array(f['PVeloc'])
	pvx = pv[0:3*Px*Py-2:3]
	#ppvx = -np.average(pvx,axis=1)
	#plt.plot(ppvx,np.arange(Py))
	plt.plot(pvx,np.arange(Np))
	plt.xlabel(r'$U_x^{p}$')
	plt.ylabel(r'depth')
	# plt.hold(True)
	
	plt.subplot(3,2,6)
	Vpx.append(pvx.sum())
	plt.plot(nn,Vpx,'*')
	plt.ylabel(r'$\sum U_x^{p}$')
	plt.xlabel(r'sum')
	# plt.hold(False)
	
	plt.tight_layout()
	name = str(num[i])+'.png'
	print(name)
	plt.savefig(name)
	plt.clf()



# plt.show()
