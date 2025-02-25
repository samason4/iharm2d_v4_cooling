import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

#to call this function (from iharm2d_v4_cooling directory):
# python files/O2_conv_cour.py ./ ./dumps7 ./dumps5 ./dumps4 ./dumps0625

# paths to .5 dumps
dumpsdir7 = sys.argv[2]
outputdir7 = sys.argv[1]
if not os.path.exists(outputdir7):
	os.makedirs(outputdir7)

# paths to .5 dumps
dumpsdir5 = sys.argv[3]
outputdir5 = sys.argv[1]
if not os.path.exists(outputdir5):
	os.makedirs(outputdir5)

# paths to .4 dumps
dumpsdir4 = sys.argv[4]
outputdir4 = sys.argv[1]
if not os.path.exists(outputdir4):
	os.makedirs(outputdir4)

# paths to .0625 dumps
dumpsdir0625 = sys.argv[5] #temp 3 instead of 4
outputdir0625 = sys.argv[1]
if not os.path.exists(outputdir0625):
	os.makedirs(outputdir0625)

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

def initial_prims(min_r, min_th, dumpsdir):
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
	ut = ucon_calc(min_r, min_th, n1, n2, ndim, dumpsdir)
	alpha = -1/(3*r**(3/2)*ut)
	u = rho**game*np.exp(kel0*(game-1))
	return [alpha, u]

def ucon_calc(min_r, min_th, n1, n2, ndim, dumpsdir):
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
	
def numerical(dumpno, min_r, min_th, dumpsdir):	
	
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
	return(u)

def analytical(dumpno, prims, dumpsdir):	
	
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
	return(u)

def find_error(dumpno, dumpsdir, want_r, want_th):
    mins = find_indices(dumpsdir, want_r, want_th)
    min_r = mins[0]
    min_th = mins[1]
    prims = initial_prims(min_r, min_th, dumpsdir)
    num = numerical(dumpno, min_r, min_th, dumpsdir)
    ana = analytical(dumpno, prims, dumpsdir)
    return(abs(num-ana))

#this part is just for the comparison line:
x = []
res = []
x2 = []
res2 = []
temp_res = 1
temp_x = 1e-6
for i in range(12):
    x.append(temp_x*temp_res**(-2))
    res.append(temp_res)
    temp_res += 2

temp_res2 = 1
temp_x2 = 1e-6
for i in range(3):
    x2.append(temp_x2*temp_res2**(-1))
    res2.append(temp_res2)
    temp_res2 += 2

#finding errors and plotting:
for i in range(201):
    cour_inv = []
    errors = []
    cour_inv.append(1/0.7)
    errors.append(find_error(i, dumpsdir7, 12, np.pi/2))
    cour_inv.append(1/0.5)
    errors.append(find_error(i, dumpsdir5, 12, np.pi/2))
    cour_inv.append(1/0.4)
    errors.append(find_error(i, dumpsdir4, 12, np.pi/2))
    #find_error(i, dumpsdir25, 12, np.pi/2)
    #find_error(i, dumpsdir125, 12, np.pi/2)
    cour_inv.append(1/0.0625)
    errors.append(find_error(i, dumpsdir0625, 12, np.pi/2))
    #find_error(i, dumpsdir03125, 12, np.pi/2)

    fig1, sub1 = plt.subplots()
    sub1.loglog(cour_inv, errors, color = 'blue', label = 'Error of Test Cooling')
    sub1.loglog(res, x, color = 'red', label = 'Line of Slope N^-2 for Comparison')
    sub1.loglog(res2, x2, color = 'g', label = 'Line of Slope N^-1 for Comparison')
    sub1.loglog(cour_inv, errors, 'bo')
    plt.xticks([], [])
    sub1.set_xticks([])
    sub1.set_xticks([], minor=True)
    sub1.set_xticks([1, 2, 4], ['2^0', '2^1','2^2'])
    plt.ylabel("Total Error")
    plt.xlabel("1 / The Courant Number")
    plt.title('Dump: {0:04d}'.format(i))
    plt.legend()
    plt.savefig(os.path.join(outputdir5,'error_vs_cour_{0:04d}.png'.format(i)))
    plt.close()
