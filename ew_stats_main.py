#from plot_norm import *
import numpy as numpy
from pylab import *
import os, sys
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
import scipy.optimize as opt

colors = get_colors()

def get_dist(z):

	'''approximate distance in pc from redshift'''

	return (z * 3e5 / 70.0) * 1e6


def emissivity_function(costheta):

	'''
	this could be modified to query an AGNSPEC spectrum, or include
	limb darkening. At the moment it just returns the argument!
	'''

	return costheta


def is_source_detected(costheta):
	'''
	apply a selection effect according to the emissivity_function
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
			#detected = True

			if (theta > max_angle):
				detected = False


		costhetas[j] = costheta

		if (theta > THRESHOLD):
			bal_flags[j] = "b"
		else:
			bal_flags[j] = "q"

	return costhetas, bal_flags



def function_to_minimize(params, ew_o_quasars, costhetas, distribution, bins):

	'''
	the function we have to minimise in order to uncover
	the intrinsic underlying 'face-on' distribution. 

	Parameters:

		params 			array-like, length 2
						mu and sigma for the gaussian 

		ew_o_quasars 	array-like 
						observed EW data for quasars 

		costhetas 		array-like
						cosines of theoretical angle distribution
	'''

	mu = params[0]
	sigma = params[1]

	#print params

	if mu <= 0 or sigma <= 0:
		return 1e50

	ewstar = distribution(mu, sigma, size=len(costhetas) )

	ew_for_test = ewstar / costhetas
	

	'''
	If the K-S statistic is small or the p-value is high, 
	then we cannot reject the hypothesis that the distributions 
	of the two samples are the same.
	'''

	# if we have different counts then we 
	normalisation = float(len(ew_o_quasars)) / float(len(costhetas))

	f_ew_for_test = normalisation * histogram(ew_for_test, bins)[0]
	f_ewo = histogram(ew_o_quasars, bins=bins)[0]

	# chi2 only really valid for counts > ~ 5, so mask others
	select = (f_ewo > 5) * (f_ew_for_test > 5)

	# is this correct?
	# df2 = f_ewo + f_ew_for_test
	df2 = f_ewo

	#chi2 = stats.chisquare(f_ewo, f_ew_for_test)

	chi2 = np.sum( (f_ewo[select] - f_ew_for_test[select])**2 / df2[select] )

	# return the reduced chi squared
	return chi2 / float(len(f_ewo[select]))



def check(ews, costhetas):

	#mu = 5.0 + (20.0 * np.random.random(size=10000))
	#sigma = (20.0 * np.random.random(size=10000))

	mu = np.arange(4,9,0.1)
	sigma = np.arange(3,20,0.1)


	chi2_min = 1e100
	params_min = None

	chi2_array = np.zeros([len(mu), len(sigma)])

	for i, m in enumerate(mu):
		for j, s in enumerate(sigma):
			chi2 = function_to_minimise([m,s], ews, costhetas)

			chi2_array[i,j] = chi2

			if chi2 < chi2_min:
				chi2_min = chi2
				params_min = (m,s)
				print i, chi2, m, s

	return chi2_array, mu, sigma


class selection:

	'''
	a class which contains a number of Boolean arrays used for easy
	selection of subsamples in the SDSS data 
	'''

	def __init__(self, data):

		redshift_lims = (3800.0 / 2800.0 - 1.0, 9200.0 / 5007.0 - 1.0)
		redshift_lims_b = (3800.0 / 1400.0 - 1.0, 9200.0 / 1700.0 - 1.0)

		self.nonbal = (data["bal_flag"] == 0) 
		self.mgbal = (data["bal_flag"] == 2)
		self.bal = (data["bal_flag"] > 0)
		self.z = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])
		self.has_o3 = (data["ew_o3"] > 0)
		self.mass = None
		self.general = self.has_o3 * self.z
		self.a = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])
		self.b = (data["z"] > redshift_lims_b[0]) * (data["z"] < redshift_lims_b[1]) * (data["ew_c4"] > 0)



class simulation:

	'''
	a class which contains a set of quantities for different angles 
	which are measures of how well that specific model did compared
	to the data. Also allows one to run the sim.
	'''

	def __init__(self, thetamin, thetamax, data, select, line_string = "ew_o3"):

		'''
		thetamin 	array-like
					array of thetamins used 

		thetamax 	array-like
					array of thetamaxes used 

		data 		dictionary
					data read in from the SDSS quasar catalog 
		'''

		shape = (len(thetamin), len(thetamax))

		self.f_bal = np.zeros(shape)		# bal frac
		self.mean = np.zeros(shape)			# mean of mock bal data
		self.std_dev = np.zeros(shape)		# std dev of mock bal data
		self.ks = np.zeros(shape)			# ks test statistic for bals
		self.ks_p_value = np.zeros(shape)	# ks p value for bals
		self.mu = np.zeros(shape)			# mean of "intrinsic" distribution
		self.sigma = np.zeros(shape)		# sigma of "intrinsic" istribution
		self.chi2 = np.zeros(shape)			# chi2 when fitting quasars
		self.thetamin = thetamin 			# angles to use
		self.thetamax = thetamax
		self.data = data
		self.line_string = line_string		# the line we are looking at
		self.bins = np.arange(0,150,2)		# bins to use

		if line_string == "ew_o3":	# use normal
			self.s_bals = select.general * select.mgbal 
			self.s_nonbals =  select.general * select.nonbal 
			self.distribution = np.random.normal

		elif line_string == "ew_c4": # use log normal
			self.s_bals = select.b * select.bal 
			self.s_nonbals = select.b * select.nonbal 
			self.distribution = np.random.lognormal


	def get_bounds(self):

		'''
		relevant bounds and guesses for different distributions
		'''

		if self.distribution == np.random.normal:
			return ((1,50),(1,50)), [5,10]
		elif self.distribution == np.random.lognormal:
			return ((1,5),(1,5)), [3,0.5]


	def run(self, select):
		'''
		run the actual simulation, using select as the booleans
		and data as the data

		populates the f_bal and ks_test arrays

		Parameters:
			select 		selection class instance
		'''

		logbins = False

		# now iterate over the thresold angles
		for i, thmin in enumerate(self.thetamin):
			for j,thmax in enumerate(self.thetamax):

				print i, j

				if thmax > thmin:

					valError = False

					#ew_obs = self.data[self.line_string][select.general*select.nonbal]
					# use sample b for the moment
					ew_obs = self.data[self.line_string][self.s_nonbals]

					NPTS =	len(ew_obs)

					# get angles above the threshold generated uniformly on sky
					# and apply to EW measurements
					costhetas_qsos, dummy = get_mock_angles(0.0, NPTS, max_angle=thmin)


					costhetas, mock_bal_flags = get_mock_angles(thmin, NPTS, max_angle=thmax)

					# do a minimization to get mu and sigma
					# for the 'intrinsic' distribution

					try:
						bounds, guess = self.get_bounds()
						minimizeObj = opt.minimize(function_to_minimize, guess,
                                               bounds = bounds, method="Powell",
					                           args = (ew_obs, costhetas_qsos, self.distribution, self.bins))
					except ValueError:
						print "Value Error!"
						valError = True

					# if either the mu or sigma returned was zero then we have a problem...
					if any(minimizeObj.x) < 0:	
						valError = True

					if valError == False:
						# copy some attributes from minimizeObj 
						self.chi2[i,j] = minimizeObj.fun
						mu, sig = minimizeObj.x

						# mock data needs to be divided by emissivity function
						mock_data = self.distribution(mu, sig, size = NPTS)
						mock_data = mock_data / emissivity_function(costhetas)

						# selections based on flags returned from angle sim
						select_mock_bals = (mock_bal_flags == "b")
						select_mock_nonbals = (mock_bal_flags == "q")

						# bal fraction- just number of objects!
						bal_frac = float(np.sum(select_mock_bals)) / float(NPTS)

						ews = self.data[self.line_string]
						if logbins: 
							mock_data = np.log10(mock_data)
							ews = np.log10(ews)

						'''
						If the K-S statistic is small or the p-value is high, 
						then we cannot reject the hypothesis that the distributions 
						of the two samples are the same.
						'''
						print i, j, thmin, thmax
						print '------------------'

						# record the values in some arrays
						self.ks_p_value[i,j] = stats.ks_2samp(mock_data[select_mock_bals], ews[self.s_bals])[1]
						self.f_bal[i,j] = bal_frac

						# store the standard dev and mean of the dataset
						self.std_dev[i,j] = np.std(mock_data[select_mock_bals])
						self.mean[i,j] = np.mean(mock_data[select_mock_bals])

						# store the attributes of the intrinsic gaussian
						self.mu[i,j] = mu
						self.sigma[i,j] = sig 


						# make a histogram with chi for this run
						ew1 = ews[self.s_nonbals]
						ew2 = self.distribution(mu, sig, size = NPTS) / costhetas_qsos

						# make some histograms 
						make_plots.individual_histogram(self.line_string, thmin, thmax, 
							                            ew1, ew2,
							                            self.chi2[i,j], self.bins)

						make_plots.bal_histogram(self.line_string, thmin, thmax, 
							                            ews[self.s_bals], mock_data[select_mock_bals],
							                            self.ks_p_value[i,j], self.f_bal[i,j], self.bins)


						print bal_frac, self.ks_p_value[i,j], self.chi2[i,j]
						print 

		return 0


	def read_from_file(self, fname="simulation.out"):

		f = open(fname, "r")

		for line in f:
			data = line.split()
			i = int(data[0])
			j = int(data[1])
			self.thetamin[i] = float(data[2])
			self.thetamax[j] = float(data[3])

			self.f_bal[i,j] = float(data[4])
			self.mean[i,j] = float(data[5])
			self.std_dev[i,j] = float(data[6])
			self.ks_p_value[i,j] = float(data[7])

		f.close()



#def function_to_minimise(A, mu, sigma, )
def write_sim_to_file(sim, fname="simulation.out"):

	f = open(fname, "w")

	f.write("# SIMULATION OUTPUT\n")
	f.write("# i, j, thmin, thmax, f_bal, mean, std_dev, ks_p_value\n")

	# create an array to write
	for i, thmin in enumerate(sim.thetamin):
		for j,thmax in enumerate(sim.thetamax):
			f.write("%i %i %.4f %.4f %8.4e %8.4e %8.4e %8.4e\n" % 
				     (i, j, thmin, thmax, sim.f_bal[i,j], sim.mean[i,j], 
			          sim.std_dev[i,j], sim.ks_p_value[i,j]) )

	f.close()



if __name__ == "__main__":

	# define angles to run sim over
	thetamin = np.arange(90,0,-5)
	thetamax = np.arange(90,0,-5)

	#thetamin = [50]
	#thetamax = [85]

	# get the data
	data = sub.get_sample("data/catalog.dat")

	# selection criteria
	# we need to be at a certain redshift and have O3 emission
	# then we select for BAL or no-BAL
	select = selection(data)
	set_pretty()

	mode = sys.argv[1]

	if mode == "sim":

		# set up sim
		sim = simulation(thetamin, thetamax, data, select, line_string = "ew_c4")

		# run sim
		sim.run(select)
		#sim.run_gauss(select)

		# make the contour plots
		make_plots.plot_contour(sim)
		#make_plots.plot_contour2(sim)

	#elif mode == "read":

	elif mode == "hist":

		make_plots.make_hist(data, select)

		mixedCase
		make_plots



