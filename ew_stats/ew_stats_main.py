#from plot_norm import *
import numpy as numpy
from pylab import *
import os, sys
from cobra_sub import smooth
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from plot_norm import *
from constants import *
from pretty import *
import scipy.stats as stats
import matplotlib.mlab as mlab
import sdss_sub_data as sub
from matplotlib.colors import LinearSegmentedColormap
cm_use=get_viridis()

def get_dist(z):

	'''approximate distance in pc from redshift'''

	return (z * 3e5 / 70.0) * 1e6


def emissivity_function(costheta):

	'''this could be modified to query an AGNSPEC spectrum'''

	return costheta


def is_source_detected(costheta):
	'''
	apply a selection effect according to the emissivity emissivity_function
	'''

	# random number
	dv = np.random.random()

	r = np.power(dv, 1.0/3.0)

	d = get_dist(r) * PARSEC

	# now we've got a random distance proportional to volume
	# get the flux assuming some luminosity (1e46?)
	flux = (1e46 * emissivity_function(costheta)) / (4.0*np.pi*d*d)

	# work out if above our detection limit which is abritrary at the moment
	d_detect = get_dist(1.0) * PARSEC
	detection_limit = 1e45 / (4.0*np.pi*d*d)

	return True


def get_mock_angles(THRESHOLD, NPTS, max_angle=None):
	'''
	generate angles according to solid angle and 
	apply selection effect 

	return flags of whether bal or not
	'''

	costhetas = np.zeros(NPTS)
	bal_flags = np.empty(NPTS, dtype="string")

	if max_angle == None:
		max_angle = 1e5

	for j in range(NPTS):

		detected = False

		while not detected:

			costheta = np.random.random()
			theta = (np.arccos(costheta) * 180.0 / np.pi)

			detected = is_source_detected(costheta)

			if (theta > max_angle):
				detected = False


		costhetas[j] = costheta

		if (theta > THRESHOLD):
			bal_flags[j] = "b"
		else:
			bal_flags[j] = "q"

	return costhetas, bal_flags


def set_subplot_ticks():
	'''
	plotting hack
	'''

	xlabel(r"$W_{\lambda}$ (\AA)", fontsize=24)
	xlim(-1,4)
	ylim(1e-3,2)
	gca().set_xticklabels(["0.1","1","10","$10^2$","$10^3$","$10^4$"])

	return 


class selection:

	'''
	a class which contains a number of Boolean arrays used for easy
	selection of subsamples in the SDSS data 
	'''

	def __init__(self, data):

		redshift_lims = (3800.0 / 2800.0 - 1.0, 9200.0 / 5007.0 - 1.0)

		self.nonbal = (data["bal_flag"] == 0) 
		self.mgbal = (data["bal_flag"] == 2)
		self.z = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])
		self.has_o3 = (data["ew_o3"] > 0)
		self.mass = None
		self.general = self.has_o3 * self.z



class simulation:

	'''
	a class which contains a set of quantities for different angles 
	which are measures of how well that specific model did compared
	to the data

	UNFINISHED
	'''

	def __init__(self, shape, NPTS):

		f_bal = np.zeros(shape)
		means = np.zeros(shape)
		std_dev = np.zeros(shape)
		ks = np.zeros(shape)
		ks_p_value = np.zeros(shape)


		






thetamins = np.arange(0,90,5)
thetamax = np.arange(45,90,5)
f_bal = np.zeros([len(thetamins), len(thetamax)])
ks_test = np.zeros([len(thetamins), len(thetamax)])



# get the data
data = sub.get_sample("../data/catalog.dat")

#some plotting modes
logbins = True
NORM  = True
labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA"]
set_pretty()
big_tick_labels(18)
colors=get_colors()


# selection criteria
# we need to be at a certain redshift and have O3 emission
# then we select for BAL or no-BAL
select = selection(data)



ews = data["ew_o3"]
if logbins: 
	ews = np.log10(ews)
	bins = np.arange(-2,4,0.1)

# now iterate over the thresold angles
for i, thmin in enumerate(thetamins):
	for j,thmax in enumerate(thetamax):

		if thmax > thmin:
			mock_data = data["ew_o3"][select.general]

			NPTS =	len(mock_data)

			# get angles above the threshold generated uniformly on sky
			# and apply to EW measurements
			costhetas, mock_bal_flags = get_mock_angles(thmin, NPTS, max_angle=thmax)

			# mock data needs to be divided by emissivity function
			mock_data = mock_data / emissivity_function(costhetas)

			# selections based on flags returned from angle sim
			select_mock_bals = (mock_bal_flags == "b")
			select_mock_nonbals = (mock_bal_flags == "q")

			# bal fraction- just number of objects!
			bal_frac = float(np.sum(select_mock_bals)) / float(NPTS)

			ews = data["ew_o3"]

			if logbins: 
				mock_data = np.log10(mock_data)
				ews = np.log10(ews)

			'''
			If the K-S statistic is small or the p-value is high, 
			then we cannot reject the hypothesis that the distributions 
			of the two samples are the same.
			'''
			print i, j
			print '------------------'

			# record the values in some arrays
			ks_test[i,j] = stats.ks_2samp(mock_data[select_mock_bals], ews[select.general*select.mgbal])[1]
			f_bal[i,j] = bal_frac
			print bal_frac, ks_test[i,j]
			print 



# below here it's just plotting...
figure(figsize=(13,7))
subplot(1,2,1)

contourf(thetamax, thetamins, np.log10(ks_test), extend='both', levels=np.arange(-20,0,0.5), cmap=cm_use)
colorbar()
contour(thetamax, thetamins, np.log10(ks_test), levels=np.log10(np.array([0.001,0.05,0.32])), c="k")
ylabel(r"$\theta_{min}$", fontsize=20)
xlabel(r"$\theta_{max}$", fontsize=20)
title('''$\log~(p$-value$)$. If the p-value is high, 
then we cannot reject the hypothesis that the distributions 
of the two samples are the same''', fontsize=14)

subplot(1,2,2)
contourf(thetamax, thetamins, f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
#contour(thetamax, thetamins, f_bal)
colorbar()
ylabel(r"$\theta_{min}$", fontsize=20)
xlabel(r"$\theta_{max}$", fontsize=20)
title("$f_{BAL}$", fontsize=20)

subplots_adjust(left=0.1,right=0.97,top=0.80)
savefig("contour.png", dpi=200)

