import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Ellipse, Circle
import matplotlib.pyplot as plt
import matplotlib.cm as cm


prefix = ""
prefix = prefix + "test_2_"
R = 10

num = np.arange(0,99+1)
nu = 0.05 




for i in range(len(num)):
	name = prefix + str(num[i]).zfill(4) + ".h5"
	print("start process ", name)
	f = h5.File(name)
	PX = np.array(f['Pposition'])
	px = PX[0:-2:3]
	py = PX[1:-1:3]
	Nx = int(f['Nx'][0])
	Ny = int(f['Ny'][0])
	vel = np.array(f['Velocity_0'])
	vx = vel[0:-2:3]
	vx = vx.reshape((Ny,Nx))
	vy = vel[1:-1:3]
	vy = vy.reshape((Ny,Nx))
	vv = np.sqrt(vx**2+vy**2)

	ax = plt.figure().add_subplot(111)
	plt.pcolor(vv,cmap='jet')
	plt.clim(vmin=0, vmax=0.05)
	# plt.colorbar()
	for j in range(np.size(PX)/6): 
		cir = Circle(xy = (px[j],py[j]), radius=R, alpha=1.0,facecolor='k')
		ax.add_patch(cir)
	plt.axis('equal')
	plt.axis('off')
	plt.tight_layout()
	name = 's'+str(i)+'.jpg'
	plt.savefig(name,dpi=500)

	plt.clf()
	plt.close()
