import numpy as np 
import matplotlib.pyplot as plt
import h5py as h5

H = 300
Ga = 27.5
D = 20
Z = np.load("Zz.npy")
zz = np.average(Z, axis=1)
temp = np.where(zz[300:]<0.1)

hf = temp[0][-1]-temp[0][0]
print("hf ", hf)
print(Ga*(hf*1.0/(D*1.0))**2)

prefix = "/media/pzhang/Elements/move-bed-tmp/macondo/"
prefix = prefix + "5e3Re_27.5Ga_0.3gap_a/"
prefix = prefix + "test_mvbed_c_"

Py = 21
Px = 160
Np = Py*Px
num = np.arange(500,924+1)
pp = []
for i in range(len(num)):
	name = prefix + str(num[i]).zfill(4) + ".h5"
	f = h5.File(name)
	Nx = np.array(f['Nx'])
	Ny = np.array(f['Ny'])
	yk = Ny[0]-hf
	p = np.array(f['Density_0'])
	p = p.reshape((Ny[0],Nx[0]))
	plt.plot(p[yk,:])
	pp.append(p[yk,:])
	name = 'p'+str(num[i])+'.png'
	plt.tight_layout()
	print(name)
	plt.savefig(name)
	plt.clf()
f = h5.File('P.h5','w')
f.create_dataset('P',data=pp)
f.close()
pp = np.array(pp)
plt.plot(np.average(pp,axis=0))
plt.show()
	