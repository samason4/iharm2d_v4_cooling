import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

#to call this function (from iharm2d_v4_cooling directory):
# python files/O2_conv_v2.py ./ ./cour=0.5_dumps ./cour=0.25_dumps ./cour=0.125_dumps ./cour=0.0625_dumps ./cour=0.03125_dumps

# paths to .5 dumps
dumpsdir5 = sys.argv[2] #second argument when you call the function (in the long email plot_density this is the "./dumps")
outputdir5 = sys.argv[1] #third argument when you call the function (in the long email plot_density this is the "./"), etc
if not os.path.exists(outputdir5):
	os.makedirs(outputdir5)

# paths to .25 dumps
dumpsdir25 = sys.argv[3]
outputdir25 = sys.argv[1]
if not os.path.exists(outputdir25):
	os.makedirs(outputdir25)

# paths to .125 dumps
dumpsdir125 = sys.argv[4]
outputdir125 = sys.argv[1]
if not os.path.exists(outputdir125):
	os.makedirs(outputdir125)

# paths to .0625 dumps
dumpsdir0625 = sys.argv[5]
outputdir0625 = sys.argv[1]
if not os.path.exists(outputdir0625):
	os.makedirs(outputdir0625)

# paths to .03125 dumps
dumpsdir03125 = sys.argv[5]
outputdir03125 = sys.argv[1]
if not os.path.exists(outputdir03125):
	os.makedirs(outputdir03125)

#function that finds the indices:
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
	
def numerical(dumpno, uarr, tarr, min_r, min_th, dumpsdir):	
	
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

def analytical(dumpno, uarr, tarr, min_r, min_th, prims, dumpsdir):	
	
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

# these are the arrays I will plot:
errors = []
resolutions = []
x = []
res = []

#finding error for the .5 run:
errors5_temp = []
for i in range(10):
	perc = i*2
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir5, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir5)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir5)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir5)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors5_temp.append(error)
error5 = 0
for i in range(10):
	error5 += abs(errors5_temp[i])
error5 = error5/11
errors.append(error5)
resolutions.append(2)

#finding error for the .25 run:
errors25_temp = []
for i in range(10):
	perc = i*2 + 20
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir25, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir25)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir25)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir25)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors25_temp.append(error)
error25 = 0
for i in range(10):
	error25 += abs(errors25_temp[i])
error25 = error25/11
errors.append(error25)
resolutions.append(4)

#finding error for the .125 run:
errors125_temp = []
for i in range(10):
	perc = i*2 + 40
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir125, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir125)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir125)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir125)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors125_temp.append(error)
error125 = 0
for i in range(10):
	error125 += abs(errors125_temp[i])
error125 = error125/11
errors.append(error125)
resolutions.append(8)

#finding error for the .0625 run:
errors0625_temp = []
for i in range(10):
	perc = i*2 + 60
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir0625, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir0625)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir0625)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir0625)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors0625_temp.append(error)
error0625 = 0
for i in range(10):
	error0625 += abs(errors0625_temp[i])
error0625 = error0625/11
errors.append(error0625)
resolutions.append(16)

#finding error for the .03125 run:
errors03125_temp = []
for i in range(10):
	perc = i*2 + 80
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir03125, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir03125)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir03125)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir03125)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors03125_temp.append(error)
error03125 = 0
for i in range(10):
	error03125 += abs(errors03125_temp[i])
error03125 = error03125/11
errors.append(error03125)
resolutions.append(32)

#this part is just for the comparison line:
temp_res = 90
temp_x = 7e-2
for i in range(27):
    temp_res += 12
    x.append(temp_x*temp_res**(-2))
    res.append(temp_res)

#actually plotting:
fig1, sub1 = plt.subplots()
sub1.loglog(resolutions, errors, color = 'blue', label = 'Error of Test Cooling')
#sub1.loglog(res, x, color = 'red', label = 'Line of Slope N^-2 for Comparison')
sub1.loglog(resolutions, errors, 'bo')
plt.xticks([], [])
sub1.set_xticks([])
sub1.set_xticks([], minor=True)
sub1.set_xticks([32, 16, 8, 4, 2], ['2^1','2^2', '2^3', '2^4', '2^5'])
plt.ylabel("Total Error")
plt.xlabel("1 / The Courant Number")
plt.title("Error vs 1/Courant Number at 96x96 resolution")
plt.legend()
plt.savefig(os.path.join(outputdir5,'error_vs_courant'))
plt.close()
