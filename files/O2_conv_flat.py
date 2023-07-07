import numpy as np
import h5py, sys, glob, psutil, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

# python files/u_time_series_flat.py ./ ./dumps9 ./dumps5 ./dumps2 ./dumps05

# paths
dumpsdir9 = sys.argv[2]
outputdir9 = sys.argv[1]
if not os.path.exists(outputdir9):
	os.makedirs(outputdir9)

dumpsdir5 = sys.argv[3]
outputdir5 = sys.argv[1]
if not os.path.exists(outputdir5):
	os.makedirs(outputdir5)

dumpsdir2 = sys.argv[4]
outputdir2 = sys.argv[1]
if not os.path.exists(outputdir2):
	os.makedirs(outputdir2)

dumpsdir05 = sys.argv[5]
outputdir05 = sys.argv[1]
if not os.path.exists(outputdir05):
	os.makedirs(outputdir05)

def initial_prims(dumpsdir, min_r, min_th):
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(0000)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
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
	
def numerical(dumpsdir, dumpno, min_r, min_th):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()
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
	return([u,t])

def analytical(dumpsdir, dumpno, prims):	
	
	# header info
	header = open(os.path.join(dumpsdir,'dump_0000{0:04d}'.format(dumpno)),'r')
	firstline = header.readline()
	header.close()
	firstline = firstline.split()

	t = float(firstline[23])
	t = '{:.3f}'.format(t)
	
	# adding u and t
	alpha = -0.2
	u0 = prims[1]
	u = u0*np.exp(alpha*float(t))
	return([u, t])

def find_error(dumpsdir, dump_min, dump_max):
	error = 0
	for i in range(dump_min, dump_max + 1):
		min_r = 5
		min_th = 5
		prims = initial_prims(dumpsdir, min_r, min_th)
		num = numerical(dumpsdir, i, min_r, min_th)[0]
		ana = analytical(dumpsdir, i, prims)[0]
		error += abs(num-ana)
	error = error/((dump_max+1)-dump_min)
	return(error)

#finding errors:
cour_inv = []
errors = []
cour_inv.append(1/0.9)
errors.append(find_error(dumpsdir9, 80, 100))
cour_inv.append(1/0.5)
errors.append(find_error(dumpsdir5, 80, 100))
cour_inv.append(1/0.2)
errors.append(find_error(dumpsdir2, 80, 100))
cour_inv.append(1/0.05)
errors.append(find_error(dumpsdir05, 80, 100))

#this is for the comparison line:
x = []
res = []
x2 = []
res2 = []
temp_res = 1
temp_x = 5e-4
for i in range(12):
    x.append(temp_x*temp_res**(-2))
    res.append(temp_res)
    temp_res += 2
temp_res2 = 1
temp_x2 = 5e-4
for i in range(12):
    x2.append(temp_x2*temp_res2**(-1))
    res2.append(temp_res2)
    temp_res2 += 2

#plotting:
fig1, sub1 = plt.subplots()
sub1.loglog(cour_inv, errors, color = 'b', label = 'Error of Test Cooling')
sub1.loglog(res, x, color = 'r', label = 'Line of Slope N^-2 for Comparison')
sub1.loglog(res2, x2, color = 'g', label = 'Line of Slope N^-1 for Comparison')
sub1.loglog(cour_inv, errors, 'bo')
plt.xticks([], [])
sub1.set_xticks([])
sub1.set_xticks([], minor=True)
sub1.set_xticks([1, 2, 4, 8, 16, 32], ['2^0', '2^1', '2^2', '2^3', '2^4', '2^5'])
plt.ylabel("Total Error")
plt.xlabel("1 / The Courant Number")
plt.title('Error vs 1/cour')
plt.legend()
plt.savefig('error_vs_cour.png')
plt.close()