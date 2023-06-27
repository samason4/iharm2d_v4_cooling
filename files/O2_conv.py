import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

#to call this function (from iharm2d_v4_cooling directory):
# python files/O2_conv.py ./ ./dumps96 ./dumps128 ./dumps256 ./dumps384

# paths to 96x96 dumps
dumpsdir96 = sys.argv[2] #second argument when you call the function (in the long email plot_density this is the "./dumps")
outputdir96 = sys.argv[1] #third argument when you call the function (in the long email plot_density this is the "./"), etc
if not os.path.exists(outputdir96):
	os.makedirs(outputdir96)

# paths to 128x128 dumps
dumpsdir128 = sys.argv[3]
outputdir128 = sys.argv[1]
if not os.path.exists(outputdir128):
	os.makedirs(outputdir128)

# paths to 256x256 dumps
dumpsdir256 = sys.argv[4]
outputdir256 = sys.argv[1]
if not os.path.exists(outputdir256):
	os.makedirs(outputdir256)

# paths to 384x384 dumps
dumpsdir384 = sys.argv[5]
outputdir384 = sys.argv[1]
if not os.path.exists(outputdir384):
	os.makedirs(outputdir384)

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

#finding errors for dump 96:
errors96_temp = []
for i in range(10):
	perc = i*2.5
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir96, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir96)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir96)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir96)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors96_temp.append(error)
error96 = 0
for i in range(10):
	error96 += abs(errors96_temp[i])
error96 = error96/11
errors.append(error96)
resolutions.append(96)

#finding errors for dump 128:
errors128_temp = []
for i in range(10):
	perc = i*2.5 + 25
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir128, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir128)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir128)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir128)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors128_temp.append(error)
error128 = 0
for i in range(10):
	error128 += abs(errors128_temp[i])
error128 = error128/11
errors.append(error128)
resolutions.append(128)

#finding errors for dump 256:
errors256_temp = []
for i in range(10):
	perc = i*2.5 + 50
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir256, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir256)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir256)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir256)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors256_temp.append(error)
error256 = 0
for i in range(10):
	error256 += abs(errors256_temp[i])
error256 = error256/11
errors.append(error256)
resolutions.append(256)

#finding errors for dump 384:
errors384_temp = []
for i in range(10):
	perc = i*2.5 + 75
	print(perc, "percent done")
	uarr_num = []
	tarr_num = []
	uarr_ana = []
	tarr_ana = []
	error = 0
	mins = find_indices(dumpsdir384, 12 + i*2, np.pi/2)
	min_r = mins[0]
	min_th = mins[1]
	prims = initial_prims(min_r, min_th, dumpsdir384)
	for i in range(181, 201): 
		numerical(i, uarr_num, tarr_num, min_r, min_th, dumpsdir384)
		analytical(i, uarr_ana, tarr_ana, min_r, min_th, prims, dumpsdir384)
	for i in range(19):
		error += abs(uarr_num[i]-uarr_ana[i])
	error = error/20
	errors384_temp.append(error)
error384 = 0
for i in range(10):
	error384 += abs(errors384_temp[i])
error384 = error384/11
errors.append(error384)
resolutions.append(384)

#this part is just for the comparison line:
temp_res = 90
temp_x = 4e-2
for i in range(27):
    temp_res += 12
    x.append(temp_x*temp_res**(-2))
    res.append(temp_res)

#actually plotting:
fig1 = plt.figure()
plt.loglog(resolutions, errors, 'bo')
plt.loglog(res, x, 'r')
plt.xlabel("resolution")
plt.ylabel("total error")
plt.title("error vs resolution")
plt.savefig(os.path.join(outputdir96,'error_vs_resolution'))
plt.close()
