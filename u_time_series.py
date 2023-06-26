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

def initial_prims(min_r, min_th):
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0000)),'r')
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
	
	#to access r:
	grid = np.loadtxt(os.path.join(dumpsdir,'grid'))
	r = grid[:,2].reshape((n1,n2))[min_r][min_th]

	# to access the other prims:
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0000)),skiprows=1)
	kel0 = prims[:,9].reshape((n1,n2))[min_r][min_th]
	rho = prims[:,0].reshape((n1,n2))[min_r][min_th]
	game = 1.333333
	m = 3
	ut = ucon_calc(min_r, min_th, n1, n2, ndim)
	alpha = -1/(3*r**(3/2)*ut)
	u = rho**game*np.exp(kel0*(game-1))
	return [alpha, u]

def ucon_calc(min_r, min_th, n1, n2, ndim):
	# reading grid file
	grid = np.loadtxt(os.path.join(dumpsdir,'grid'))
	gcov = grid[:,24:].reshape((n1, n2, ndim, ndim))
	lapse = grid[:,7].reshape((n1,n2))

	# loading prims
	prims1 = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0000)),skiprows=1)
	u1 = prims1[:,2].reshape(n1,n1)[min_r][min_th]
	u2 = prims1[:,3].reshape(n1,n1)[min_r][min_th]
	u3 = prims1[:,4].reshape(n1,n1)[min_r][min_th]

	# calculating ut
	qsq = gcov[min_r][min_th][1][1]*u1*u1+gcov[min_r][min_th][2][2]*u2*u2+gcov[min_r][min_th][3][3]*u3*u3+2*(gcov[min_r][min_th][1][2]*u1*u2+gcov[min_r][min_th][1][3]*u1*u3+gcov[min_r][min_th][2][3]*u2*u3)
	gamma = (1+qsq)**(0.5)
	ut = gamma/lapse[min_r][min_th]
	return ut
	
def numerical(dumpno, uarr, tarr, min_r, min_th):	
	
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
	
	# loading prims
	prims = np.loadtxt(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),skiprows=1)
	kel0 = prims[:,9].reshape((n1,n2))
	rho = prims[:,0].reshape((n1,n2))
	game = 1.333333
	u = rho[min_r][min_th]**game*np.exp(kel0[min_r][min_th]*(game-1))
	uarr.append(u)
	tarr.append(t)

def analytical(dumpno, uarr, tarr, min_r, min_th, prims):	
	
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
	
	# adding u and t
	alpha = prims[0]
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
mins = find_indices()
min_r = mins[0]
min_th = mins[1]
prims = initial_prims(min_r, min_th)
print("finding numerical and analytical solutions...")
for i in range(201): 
	numerical(i, uarr_num, tarr_num, min_r, min_th)
	analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims)

fig1 = plt.figure()
plt.plot(tarr_ana, uarr_ana, 'b')
plt.plot(tarr_num, uarr_num, 'r.')
plt.xticks(range(0, 201, 50))
plt.xlabel("time")
plt.ylabel("electron internal energy")
plt.title("numerical->red, analytical->blue")
plt.savefig(os.path.join(outputdir,'u_vs_t'))
plt.close()

print("finding the error...")
for i in range(201):
	diff.append(abs(uarr_num[i]-uarr_ana[i]))
fig2 = plt.figure()
plt.plot(tarr_ana, diff, 'b')
plt.xticks(range(0, 201, 50))
plt.xlabel("time")
plt.ylabel("error")
plt.title("total error in internal energy over time")
plt.savefig(os.path.join(outputdir,'error_u_vs_t'))
plt.close()