import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

#to call this function (from iharm2d_v4_cooling directory):
# python files/O2_conv_res_and_cour.py ./ ./dumps9 ./dumps5 ./dumps2 ./dumps05 ./dumps96 ./dumps128 ./dumps256 ./dumps384

# paths to .9 dumps
dumpsdir9 = sys.argv[2]
outputdir9 = sys.argv[1]
if not os.path.exists(outputdir9):
	os.makedirs(outputdir9)

# paths to .5 dumps
dumpsdir5 = sys.argv[3]
outputdir5 = sys.argv[1]
if not os.path.exists(outputdir5):
	os.makedirs(outputdir5)

# paths to .2 dumps
dumpsdir2 = sys.argv[4]
outputdir2 = sys.argv[1]
if not os.path.exists(outputdir2):
	os.makedirs(outputdir2)

# paths to .05 dumps
dumpsdir05 = sys.argv[5]
outputdir05 = sys.argv[1]
if not os.path.exists(outputdir05):
	os.makedirs(outputdir05)

# paths to 96 dumps
dumpsdir96 = sys.argv[6]
outputdir96 = sys.argv[1]
if not os.path.exists(outputdir96):
	os.makedirs(outputdir96)

# paths to 128 dumps
dumpsdir128 = sys.argv[7]
outputdir128 = sys.argv[1]
if not os.path.exists(outputdir128):
	os.makedirs(outputdir128)

# paths to 256 dumps
dumpsdir256 = sys.argv[8]
outputdir256 = sys.argv[1]
if not os.path.exists(outputdir256):
	os.makedirs(outputdir256)
	
# paths to 384 dumps
dumpsdir384 = sys.argv[9]
outputdir384 = sys.argv[1]
if not os.path.exists(outputdir384):
	os.makedirs(outputdir384)

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
	dt = firstline[38]

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
	return([u, dt])

def analytical(dumpno, prims, dumpsdir):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
	metric = firstline[9]
	dt = firstline[38]

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
	return([u, dt])

def find_error(dumpsdir, want_r, want_th, errors, dts):
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	mins = find_indices(dumpsdir, want_r, want_th)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir)
	num = numerical(199, min_r, min_th, dumpsdir)
	ana = analytical(199, prims, dumpsdir)
	error = abs(ana[0]-num[0])
	errors.append(error)
	dts.append(float(num[1]))
	#print(float(num[1]))
    
# these are the arrays I will plot:
errors_cour = []
dts_cour = []
errors_res = []
dts_res = []
x = []
res = []
x2 = []
res2 = []

#finding errors:
find_error(dumpsdir9, 12, np.pi/2, errors_cour, dts_cour)
find_error(dumpsdir5, 12, np.pi/2, errors_cour, dts_cour)
find_error(dumpsdir2, 12, np.pi/2, errors_cour, dts_cour)
find_error(dumpsdir05, 12, np.pi/2, errors_cour, dts_cour)
find_error(dumpsdir96, 12, np.pi/2, errors_res, dts_res)
find_error(dumpsdir128, 12, np.pi/2, errors_res, dts_res)
find_error(dumpsdir256, 12, np.pi/2, errors_res, dts_res)
find_error(dumpsdir384, 12, np.pi/2, errors_res, dts_res)

#this part is just for the comparison line:
"""temp_res = 1
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
    temp_res2 += 2"""

#actually plotting:
fig1, sub1 = plt.subplots()
print('dts_res: ', dts_res)
print('errors_res: ', errors_res)
print('dts_cour: ', dts_cour)
print('errors_cour: ', errors_cour)
sub1.loglog(dts_res, errors_res, color = 'b', label = 'Error from varying resolution')
sub1.loglog(dts_cour, errors_cour, color = 'y', label = 'Error from varying courant')
#sub1.loglog(res, x, color = 'r', label = 'Line of Slope N^-2 for Comparison')
#sub1.loglog(res2, x2, color = 'g', label = 'Line of Slope N^-1 for Comparison')
sub1.loglog(dts_res, errors_res, 'bo')
sub1.loglog(dts_cour, errors_cour, 'yo')
"""plt.xticks([], [])
sub1.set_xticks([])
sub1.set_xticks([], minor=True)
sub1.set_xticks([0.05, 0.002], ['0.05', '0.002'])"""
plt.ylabel("Total Error")
plt.xlabel("dt")
plt.title("Error vs dt with different courant numbers and resolutions")
plt.legend()
plt.savefig(os.path.join(outputdir5,'error_vs_dt'))
plt.close()
