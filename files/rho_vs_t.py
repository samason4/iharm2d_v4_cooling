import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

#to call this function (from iharm2d_v4_cooling directory):
# python files/rho_vs_t.py ./ ./dumps

# paths to .5 dumps
dumpsdir9 = sys.argv[2]
outputdir9 = sys.argv[1]
if not os.path.exists(outputdir9):
	os.makedirs(outputdir9)

def find_indices(dumpsdir, r_want, th_want):

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
	min_r  = np.argmin(np.fabs(r[:,0] - r_want))
	min_th = np.argmin(np.fabs(th[min_r,:] - th_want))
	return [min_r, min_th]
	
def log_rho(dumpno, min_r, min_th, dumpsdir):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	n1 = int(firstline[11])
	n2 = int(firstline[12])
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	rho = prims[:,0].reshape((n1,n2))[min_r][min_th]
	print(rho)
	return(-np.log10(rho))

def find_rho(dumpsdir, want_r, want_th, nums, ts):
    small = 10
    large = 0.00000000
    mins = find_indices(dumpsdir, want_r, want_th)
    min_r = mins[0]
    min_th = mins[1]
    for i in range(201):
    	nume = log_rho(i, min_r, min_th, dumpsdir)
	    nums.append(nume)
        ts.append(5*i)
        if(small>nume):
            small = nume
    	if(large<nume):
            large = nume
    print(large-small)



#finding errors:
ts = []
errors = []
find_rho(dumpsdir9, 22, np.pi/2, errors, ts)

#plotting:
fig1, sub1 = plt.subplots()
sub1.plot(ts, errors, color = 'b')
plt.ylabel("-Log(rho) at r=22 and theta=pi/2")
plt.xlabel("time")
plt.title('-Log(rho) vs time')
plt.savefig('rho_vs_time.png')
plt.close()
