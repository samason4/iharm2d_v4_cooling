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

# function to parallelize plotting
def run_parallel(function, dlist,	nthreads):
	pool = mp.Pool(nthreads)
	pool.map_async(function, dlist).get(720000)
	pool.close()
	pool.join()

# function to generate 4-potential from magnetic field
def plotting_bfield_lines(ax, B1, B2, xp, zp, n1, n2, dx1, dx2, gdet, nlines=20):
	AJ_phi = np.zeros([n1,n2]) 
	for j in range(n2):
			for i in range(n1):
					AJ_phi[i,j] = (np.trapz(gdet[:i,j]*B2[:i,j],dx=dx1) - np.trapz(gdet[i,:j]*B1[i,:j],dx=dx2))
	AJ_phi -=AJ_phi.min()
	levels = np.linspace(0,AJ_phi.max(),nlines*2)
	ax.contour(xp, zp, AJ_phi, levels=levels, colors='k')

def plotting(dumpno):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()

	madtype = int(firstline[0])
	rin = float(firstline[2])	
	rmax = float(firstline[3])
	n1 = int(firstline[11])
	n2 = int(firstline[12])

	t = '{:.3f}'.format(t)
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	u = prims[:,1].reshape((n1,n2))
	max = np.argmax(u, axis=0, out = none)
	index1 = 0
	for i in range(n2):
		if(max > n1):
			index1++
			max -= n1
	u_max = u[index1][max]
	print(umax)
	
