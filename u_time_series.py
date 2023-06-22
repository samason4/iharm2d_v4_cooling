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

def plotting(dumpno):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
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
	for i in range(n1):
		for j in range(n2):
			r[i][j] -= 12
			r[i][j] = abs(r.all())
			th[i][j] -= np.pi/2
			th[i][j] = abs(th.all())
	minarr_r = np.argmin(r, axis=0)
	min_r = minarr_r[0]
	minarr_th = np.argmin(th, axis=1)
	min_th = minarr_th[min_r]
	print("coords:", (min_r, min_th))
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	u = prims[:,1].reshape((n1,n2))

#actually call the function:
plotting(0000)
