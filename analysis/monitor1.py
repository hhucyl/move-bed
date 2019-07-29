import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

prefix = "/media/pzhang/Elements/move-bed-tmp/macondo/"
prefix = prefix + "0.01v_0.01nu_40l_30x_10y/1/"
prefix = prefix + "test_mvbed_s_"

num = np.arange(int(sys.argv[1])+1)
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
	
	plt.subplot(3,1,1)
	nn.append(num[i])
	vx = np.array(vx)
	idx = int(Nx[0]/2)
	vvvx = vx[:,idx]
	Vx.append(np.abs(vvvx.sum()))
	plt.plot(nn,Vx,'*')
	plt.ylabel(r'$\sum u_x$')
	plt.xlabel(r'num')
	# plt.hold(False)


	plt.subplot(3,1,2)
	pv = np.array(f['PVeloc'])
	Np = int(int(pv.size)/6)
	pvx = pv[0:3*Np-2:3]
	#ppvx = -np.average(pvx,axis=1)
	#plt.plot(ppvx,np.arange(Py))
	plt.plot(pvx,np.arange(Np))
	plt.xlabel(r'$U_x^{p}$')
	plt.ylabel(r'depth')
	# plt.hold(True)
	
	plt.subplot(3,1,3)
	Vpx.append(pvx.sum())
	plt.plot(nn,Vpx,'*')
	plt.ylabel(r'$\sum U_x^{p}$')
	plt.xlabel(r'sum')
	# plt.hold(False)
	
	plt.tight_layout()
	name = 's'+str(num[i])+'.png'
	print(name)
	plt.savefig(name)
	plt.clf()



# plt.show()
