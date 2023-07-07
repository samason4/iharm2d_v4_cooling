import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

# python files/u_time_series_flat.py ./dumps ./

# paths
dumpsdir = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir):
	os.makedirs(outputdir)

def initial_prims(min_r, min_th):
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0000)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	metric = firstline[4]
	n1 = int(firstline[6])
	n2 = int(firstline[7])

	t = float(firstline[23])
	t = '{:.3f}'.format(t)

	# to access the other prims:
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0000)),skiprows=1)
	kel0 = prims[:,9].reshape((n1,n2))[min_r][min_th]
	rho = prims[:,0].reshape((n1,n2))[min_r][min_th]
	game = 1.333333
	alpha = -0.2
	u = rho**game*np.exp(kel0*(game-1))
	return [alpha, u]
	
def numerical(dumpno, uarr, tarr, min_r, min_th):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	metric = firstline[4]
	n1 = int(firstline[6])
	n2 = int(firstline[7])

	t = float(firstline[23])
	t = '{:.3f}'.format(t)
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	kel0 = prims[:,9].reshape((n1,n2))
	rho = prims[:,0].reshape((n1,n2))
	game = 1.333333
	u = rho[min_r][min_th]**game*np.exp(kel0[min_r][min_th]*(game-1))
	uarr.append(u)
	tarr.append(t)

def analytical(dumpno, uarr, tarr, prims):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	metric = firstline[4]

	t = float(firstline[23])
	t = '{:.3f}'.format(t)
	
	# adding u and t
	alpha = -0.2
	u0 = prims[1]
	u = u0*np.exp(alpha*float(t))
	uarr.append(u)
	tarr.append(t)

#actual plotting:
uarr_num = []
tarr_num = []
uarr_ana = []
tarr_ana = []
diff = []
min_r = 5
min_th = 5
prims = initial_prims(min_r, min_th)
print("finding numerical and analytical solutions...")
for i in range(101): 
	numerical(i, uarr_num, tarr_num, min_r, min_th)
	analytical(i, uarr_ana, tarr_ana, prims)

fig1 = plt.figure()
plt.plot(tarr_ana, uarr_ana, 'bo')
plt.plot(tarr_num, uarr_num, 'r.')
plt.xticks(range(0, 101, 25))
plt.xlabel("time")
plt.ylabel("electron internal energy")
plt.title("numerical->red, analytical->blue")
plt.savefig(os.path.join(outputdir,'u_vs_t'))
plt.close()

print("finding the error...")
for i in range(101):
	diff.append(abs(uarr_num[i]-uarr_ana[i]))
fig2 = plt.figure()
plt.plot(tarr_ana, diff, 'b')
plt.xticks(range(0, 101, 50))
plt.xlabel("time")
plt.ylabel("error")
plt.title("total error in internal energy over time")
plt.savefig(os.path.join(outputdir,'error_u_vs_t'))
plt.close()