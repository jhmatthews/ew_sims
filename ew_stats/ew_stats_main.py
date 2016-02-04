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
from scipy.optimize import curve_fit
import make_plots

colors = get_colors()

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

# Define model function to be used to fit to the data above:
def gauss(x, *p):
	
	'''
	gaussian function 
	'''

	A, mu, sigma = p

	return A * np.exp(-(x-mu)**2/(2.*sigma**2))


def fit_histogram(ews, max_angle):

	'''
	get fit parameters for EW for the quasars. This is our best
	guess at an intrinsic distribution.
	'''

	# first get random angles just for some quasars
	NPTS = len(ews)

	if max_angle > 0.0:
		costhetas, flags = get_mock_angles(0.0, NPTS, max_angle=max_angle)
	else:
		costhetas = np.ones(NPTS)

	# intrinsic distribution will be corrected by costhetas
	intrinsic = costhetas * ews
	log_intrinsic = np.log10(intrinsic)

	# make a histogram of log EW
	n, bins = np.histogram(log_intrinsic, normed=True, bins=np.arange(-1,3,0.05))
	
	# estimate and error, root N stats
	n_error = np.sqrt(n)

	# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
	p0 = [1., 0., 1.]

	# fit with gaussian
	bin_centres = 0.5*(bins[1:] + bins[:-1])
	coeff, var_matrix = curve_fit(gauss, bin_centres, n, p0=p0)

	return coeff, intrinsic




def set_subplot_ticks():
	'''
	plotting hack
	'''

	xlabel(r"$W_{\lambda}$ (\AA)", fontsize=24)
	xlim(-1,4)
	ylim(1e-3,2)
	gca().set_xticklabels(["0.1","1","10","$10^2$","$10^3$","$10^4$"])

	return 


def mock_data_from_gauss(NPTS, *coeff):

	A, loc, scale = coeff 

	return np.random.normal(loc=loc, scale=scale, size=NPTS)









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
	to the data. Also allows one to run the sim.
	'''

	def __init__(self, thetamin, thetamax, data):

		'''
		thetamin 	array-like
					array of thetamins used 

		thetamax 	array-like
					array of thetamaxes used 

		data 		dictionary
					data read in from the SDSS quasar catalog 
		'''

		shape = (len(thetamin), len(thetamax))

		self.f_bal = np.zeros(shape)
		self.mean = np.zeros(shape)
		self.std_dev = np.zeros(shape)
		self.ks = np.zeros(shape)
		self.ks_p_value = np.zeros(shape)
		self.thetamin = thetamin
		self.thetamax = thetamax

		self.data = data

	def run(self, select):
		'''
		run the actual simulation, using select as the booleans
		and data as the data

		populates the f_bal and ks_test arrays
		'''
		logbins = False

		# now iterate over the thresold angles
		for i, thmin in enumerate(self.thetamin):
			for j,thmax in enumerate(self.thetamax):

				if thmax > thmin:
					mock_data = self.data["ew_o3"][select.general]

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

					ews = self.data["ew_o3"]

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
					self.ks_p_value[i,j] = stats.ks_2samp(mock_data[select_mock_bals], ews[select.general*select.mgbal])[1]
					self.f_bal[i,j] = bal_frac

					# store the standard dev and mean of the dataset
					self.std_dev[i,j] = np.std(mock_data[select_mock_bals])
					self.mean[i,j] = np.mean(mock_data[select_mock_bals])

					figure()
					hist(np.log10(ews[select.general*select.mgbal]), normed=True, facecolor=colors[0], alpha=0.5, bins=np.arange(-2,3,0.1))
					hist(np.log10(mock_data[select_mock_bals]), normed=True, facecolor=colors[1], alpha=0.5, bins=np.arange(-2,3,0.1))
					savefig("hist_%i_%i.png" % (thmin, thmax), dpi=100)
					clf()

					print bal_frac, self.ks_p_value[i,j]
					print 

		return 0


	#def function_to_minimise(A, mu, sigma, )



	def run_gauss(self, select):
		'''
		run the actual simulation, using select as the booleans
		and data as the data

		populates the f_bal and ks_test arrays
		'''
		logbins = False

		# now iterate over the thresold angles
		for i, thmin in enumerate(self.thetamin):
			for j,thmax in enumerate(self.thetamax):

				print i, j

				if thmax > thmin:
					
					ews = self.data["ew_o3"][select.general]
					ew_bals = self.data["ew_o3"][select.general*select.mgbal]

					# fit the histogram yeh
					print "fitting"
					coeff, intrinsic = fit_histogram(ews, thmin)
					print "successful fit"

					NPTS =	len(ews)

					# get angles above the threshold generated uniformly on sky
					# and apply to EW measurements
					costhetas, mock_bal_flags = get_mock_angles(thmin, NPTS, max_angle=thmax)

					# mock data needs to be divided by emissivity function
					# mock_data = self.data["ew_o3"][select.general]
					mock_data = 10.0**mock_data_from_gauss(NPTS, *coeff)
					mock_data = mock_data / emissivity_function(costhetas)

					# selections based on flags returned from angle sim
					select_mock_bals = (mock_bal_flags == "b")
					select_mock_nonbals = (mock_bal_flags == "q")

					# bal fraction- just number of objects!
					bal_frac = float(np.sum(select_mock_bals)) / float(NPTS)


					'''
					If the K-S statistic is small or the p-value is high, 
					then we cannot reject the hypothesis that the distributions 
					of the two samples are the same.
					'''
					print i, j
					print '------------------'

					# record the values in some arrays
					self.ks_p_value[i,j] = stats.ks_2samp(mock_data[select_mock_bals], ew_bals)[1]
					self.f_bal[i,j] = bal_frac

					# store the standard dev and mean of the dataset
					self.std_dev[i,j] = np.std(mock_data[select_mock_bals])
					self.mean[i,j] = np.mean(mock_data[select_mock_bals])

					figure()
					bins = np.arange(-2,3,0.1)
					ALPHA = 0.3
					hist(np.log10(ew_bals), normed=True, facecolor=colors[0], alpha=ALPHA, bins=np.arange(-2,3,0.1), label="BAL")
					hist(np.log10(mock_data[select_mock_bals]), normed=True, facecolor=colors[1], alpha=ALPHA, bins=np.arange(-2,3,0.1), label="Mock BAL")
					hist(np.log10(intrinsic), normed=True, facecolor=colors[3], alpha=ALPHA, bins=np.arange(-2,3,0.1), label="Intrinsic")
					hist(np.log10(ews), normed=True, facecolor=colors[4], alpha=ALPHA, bins=np.arange(-2,3,0.1), label="BAL")
					binc = 0.5*(bins[1:] + bins[:-1])
					plot(binc, gauss(binc, *coeff), c="r", linewidth=2, label="Intrinsic fit")
					float_legend()
					savefig("hist_%i_%i.png" % (thmin, thmax), dpi=100)
					clf()

					print bal_frac, self.ks_p_value[i,j]
					print 

		return 0




if __name__ == "__main__":

	# define angles to run sim over
	thetamin = np.arange(0,90,5)
	thetamax = np.arange(0,90,5)

	# get the data
	data = sub.get_sample("../data/catalog.dat")

	# selection criteria
	# we need to be at a certain redshift and have O3 emission
	# then we select for BAL or no-BAL
	select = selection(data)

	# set up sim
	sim = simulation(thetamin, thetamax, data)

	# run sim
	#sim.run(select)
	sim.run_gauss(select)

	# make the contour plots
	make_plots.plot_contour(sim)
	make_plots.plot_contour2(sim)
