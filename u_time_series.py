import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

# paths
dumpsdir = sys.argv[1]
outputdir = sys.argv[2]
if not os.path.exists(outputdir):
	os.makedirs(outputdir)

def plotting(dumpno, uarr, tarr, min_r, min_th):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	metric = firstline[9]
	n1 = int(firstline[11])
	n2 = int(firstline[12])

	# if electron heating was enabled
	if len(firstline) > 38:
		ndim = int(firstline[27])
		if metric == 'FMKS':
			t = float(firstline[37])
		elif metric == 'MKS':
			t = float(firstline[34])
	# if electron heating was not enabled
	else:
		ndim = int(firstline[22])	
		if metric == 'FMKS':
			t = float(firstline[32])
		elif metric == 'MKS':
			t = float(firstline[29])

	t = '{:.3f}'.format(t)

	print(dumpno)
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	kel0 = prims[:,9].reshape((n1,n2))
	rho = prims[:,0].reshape((n1,n2))
	game = 1.333333
	u = rho[min_r][min_th]**game*np.exp(kel0[min_r][min_th]*(game-1))
	uarr.append(u)
	tarr.append(t)

#function that finds the indices:
def find_indices():
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	metric = firstline[9]
	n1 = int(firstline[11])
	n2 = int(firstline[12])

	# if electron heating was enabled
	if len(firstline) > 38:
		ndim = int(firstline[27])
		if metric == 'FMKS':
			t = float(firstline[37])
		elif metric == 'MKS':
			t = float(firstline[34])
	# if electron heating was not enabled
	else:
		ndim = int(firstline[22])	
		if metric == 'FMKS':
			t = float(firstline[32])
		elif metric == 'MKS':
			t = float(firstline[29])

	t = '{:.3f}'.format(t)

	#load grid stuff
	grid = np.loadtxt(os.path.join(dumpsdir,'grid'))
	r = grid[:,2].reshape((n1,n2))
	th = grid[:,3].reshape((n1,n2))

	#find the right indices
	min_r  = np.argmin(np.fabs(r[:,0] - 12))
	min_th = np.argmin(np.fabs(th[min_r,:] - np.pi/2))
	return [min_r, min_th]

#actual plotting:
uarr = []
tarr = []
mins = find_indices()
min_r = mins[0]
min_th = mins[1]
print(mins)
for i in range(201): 
	plotting(i, uarr, tarr, min_r, min_th)
plt.plot(tarr, uarr, 'b')
plt.xticks(range(0, 201, 50))
plt.xlabel("time")
plt.ylabel("electron internal energy")
plt.savefig(os.path.join(outputdir,'internal_energy_vs_time'))
plt.close()
