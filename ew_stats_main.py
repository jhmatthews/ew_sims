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


def is_source_detected(flux):
	'''
	apply a selection effect according to the emissivity_function
	'''

	# random number
	detection_limit = 10.0 ** (-12.3)

	detected = (flux > detection_limit)

	return detected


def get_mock_angles(THRESHOLD, NPTS, fluxes, max_angle=None):
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

			# select a random luminosity
			flux = fluxes[np.random.randint(0,len(fluxes))] * costheta

			detected = is_source_detected(flux)
			#detected = True

			if (theta > max_angle):
				detected = False

			#detected = True


		costhetas[j] = costheta

		if (theta > THRESHOLD):
			bal_flags[j] = "b"
		else:
			bal_flags[j] = "q"

	return costhetas, bal_flags


def get_qso_angles(NPTS, thmin, thmax):
	'''
	generate angles according to solid angle and 
	apply selection effect 

	return flags of whether bal or not
	'''

	costhetas = np.zeros(NPTS)


	for j in range(NPTS):

		detected = False

		while not detected:

			costheta = np.random.random()
			theta = (np.arccos(costheta) * 180.0 / np.pi)

			#detected = is_source_detected(costheta)
			detected = True

			if (theta > thmin) and (theta < thmax):
				detected = False


		costhetas[j] = costheta

	return costhetas

def get_bal_flags(angles, thmin, thmax):

	bal_flags = np.empty(len(angles), dtype="string")

	for i, a in enumerate(angles):

		theta = (np.arccos(a) * 180.0 / np.pi) 

		if theta > thmin and theta < thmax:
			bal_flags[i] = "b"
		else:
			bal_flags[i] = "q"

	return bal_flags

def function_to_minimize(params, ew_o_quasars, costhetas, distribution, bins):

	'''
	The reduced chi2 function we have to minimise in order to uncover
	the intrinsic underlying 'face-on' distribution.

	Parameters:

		params 			array-like, length 2
						mu and sigma for the gaussian 

		ew_o_quasars 	array-like 
						observed EW data for quasars 

		costhetas 		array-like
						cosines of theoretical angle distribution

		distribution 	function
						shape of intrinsic distribution

		bins 			array-like 
						array of bins 

	Returns:
		chi2/dof 		float 
						reduced chi2 for this model
	'''

	# get our model parameters 
	mu = params[0]
	sigma = params[1]

	# only allow positive values of mu and sigma
	if mu <= 0 or sigma <= 0:
		return 1e50

	# intrinsic distribution
	ewstar = distribution(mu, sigma, size=len(costhetas) )

	# model distribution
	ew_for_test = ewstar / costhetas

	# if we have different counts then we need to renormalise
	normalisation = float(len(ew_o_quasars)) / float(len(costhetas))
	#normalisation = 1.0

	# apply the normalisation 
	f_ew_for_test = normalisation * histogram(ew_for_test, bins=bins)[0]
	f_ewo = histogram(ew_o_quasars, bins=bins)[0]

	# chi2 only really valid for counts > ~ 5, so mask others
	select = (f_ewo > 0) * (f_ew_for_test > 0)
	#select = (f_ewo > -1)

	# the variance is sigma**2, and sigma = sqrt(N). N=f_ewo, so variance is:
	df2 = f_ewo

	# now we calculate chi squared
	chi2_arr = (f_ewo[select] - f_ew_for_test[select])**2
	chi2_arr /= df2[select]
	chi2 = np.sum(chi2_arr)

	# degrees of freedom
	dof = len(f_ewo[select]) - 3


	#print chi2, dof

	# return the reduced chi squared
	return chi2/dof



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

	def __init__(self, data, hst_map=None):

		redshift_lims = (3800.0 / 2800.0 - 1.0, 9200.0 / 5007.0 - 1.0)
		redshift_lims_b = (3800.0 / 1550.0 - 1.0, 9200.0 / 2800.0 - 1.0)

		self.nonbal = (data["bal_flag"] == 0) 
		self.mgbal = (data["bal_flag"] == 2)
		
		if hst_map == None:
			self.bal = (data["bal_flag"] > 0)
		else:
			self.bal = np.zeros(len(data["bal_flag"]))
			self.bal[hst_map] = 1
			self.bal = self.bal.astype(bool)

		self.core = (data["radio_flag"] == 1)
		self.lobe = (data["radio_flag"] == 2)
		self.z = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])
		self.has_o3 = (data["ew_o3"] > 0)
		self.mass = None
		self.general = self.has_o3 * self.z
		self.a = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1]) * self.has_o3
		self.b = (data["z"] > redshift_lims_b[0]) * (data["z"] < redshift_lims_b[1]) * (data["ew_c4"] > 0)

		dd = get_dist(data["z"])

		abs_mag = data["app"] + 5 - (5.0 * np.log10(dd))
		self.R11 = (data["z"] > 0.01) * (data["z"] < 0.8) * (data["app"] < 19.1) * (abs_mag < 22.1) * (data["SNR"] > 5)


class simulation:

	'''
	a class which contains a set of quantities for different angles 
	which are measures of how well that specific model did compared
	to the data. Also allows one to run the sim.
	'''

	def __init__(self, thetamin, thetamax, data, select, line_string = "ew_o3", mode="max", source="sdss", saves=None):

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
		self.mean_qsos = np.zeros(shape)	# mean of mock nonbal data
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
		self.mode = mode
		self.source = source
		self.saves = saves

		if line_string == "ew_o3":	# use normal
			if source == "sdss":
				self.s_bals = select.general * select.mgbal 
			elif source == "hst":
				self.s_bals = select.general * select.bal

			self.s_nonbals =  select.general * select.nonbal 
			self.distribution = np.random.lognormal
			self.R11 = select.R11


		elif line_string == "ew_c4": # use log normal
			self.s_bals = select.b * select.bal 
			self.s_nonbals = select.b * select.nonbal 
			self.distribution = np.random.lognormal

		os.system("mkdir examples_%s_%s_%s" % (self.line_string, self.mode, self.source) )


	def get_bounds(self):

		'''
		relevant bounds and guesses for different distributions
		'''

		if self.distribution == np.random.normal:
			return ((1,50),(1,50)), [5,10]
		elif self.distribution == np.random.lognormal:
			return ((-1,50),(-1,50)), [1,1]


	def run(self, select):
		'''
		run the actual simulation, using select as the booleans
		and data as the data

		populates the f_bal and ks_test arrays

		Parameters:
			select 		selection class instance
		'''

		logbins = False

		ew_obs = self.data[self.line_string][self.s_nonbals]
		fluxes = self.data["L5100"][self.s_nonbals]

		NPTS =	len(ew_obs)

		if self.mode == "nomax":
			costhetas, dummy = get_mock_angles(0.0, NPTS, fluxes, max_angle=90.0)

		# place to write the simulation output
		output_file = open("FO_simulation_%s_%s_%s.out" % (self.line_string, self.mode, self.source), "w")

		# now iterate over the thresold angles
		for i, thmin in enumerate(self.thetamin):
			for j,thmax in enumerate(self.thetamax):

				print i, j

				diff = thmax - thmin

				if diff >= 4.0:		# minimum opening angle of 4 degrees

					valError = False

					if self.mode == "max":
					# get angles above the threshold generated uniformly on sky
					# and apply to EW measurements
						costhetas_qsos, dummy = get_mock_angles(0.0, NPTS, fluxes, max_angle=thmin)
						costhetas, mock_bal_flags = get_mock_angles(thmin, NPTS, fluxes, max_angle=thmax)

					elif self.mode == "faceon":
						costhetas, mock_bal_flags = get_mock_angles(thmin, NPTS, fluxes, max_angle=thmax)

					elif self.mode == "nomax":
						costhetas_qsos = get_qso_angles(NPTS, thmin, thmax)
						mock_bal_flags = get_bal_flags(costhetas, thmin, thmax)


					# do a minimization to get mu and sigma
					# for the 'intrinsic' distribution
					# only do this if the mode is not "faceon"
					if self.mode != "faceon":
						try:
							bounds, guess = self.get_bounds()
							minimizeObj = opt.minimize(function_to_minimize, guess,
                                               bounds = bounds, method="Powell",
					                           args = (ew_obs, costhetas_qsos, self.distribution, self.bins))
						except ValueError:
							print "Value Error!"
							valError = True

						# if either the mu or sigma returned was zero then we have a problem...
						if any(minimizeObj.x) <= 0:	
							valError = True



					if valError == False:

						# get the data that we use
						#ews = self.data[self.line_string]
						ews = ew_obs

						if self.mode != "faceon":
							# copy some attributes from minimizeObj 
							self.chi2[i,j] = minimizeObj.fun
							mu, sig = minimizeObj.x

							# mock data needs to be divided by emissivity function
							mock_data = self.distribution(mu, sig, size = NPTS)
						
						else:
							mock_data = ews
							self.chi2[i,j] = -999		 # dummy value
							mu, sig = -999,-999

						print len(mock_data), NPTS

						# mock data needs to be divided by emissivity function
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
						self.mean_qsos[i,j] =  np.mean(mock_data[select_mock_nonbals])

						# store the attributes of the intrinsic gaussian
						self.mu[i,j] = mu
						self.sigma[i,j] = sig 


						#make a histogram with chi for this run
						#if thmin == 84 and thmax == 90:
						#	ew1 = ews[self.s_nonbals]
						#	ew2 = self.distribution(mu, sig, size = NPTS) / costhetas_qsos

						#	make_plots.individual_scatter(self, thmin, thmax, ew1, ew2, self.chi2[i,j], self.bins)
						#	sys.exit()

						if self.saves != None:
							for k in range(len(self.saves)):
								if thmin == self.saves[k][0] and thmax == self.saves[k][1]:
									# make some histograms 
									make_plots.individual_histogram(self, thmin, thmax, 
									                            ew1, ew2,
									                            self.chi2[i,j], self.bins)

									make_plots.bal_histogram(self, thmin, thmax, 
									                            ews[self.s_bals], mock_data[select_mock_bals],
									                            self.mean[i,j], self.f_bal[i,j], self.bins)


						print bal_frac, self.ks_p_value[i,j], self.chi2[i,j]

						#output_file.write("%i %i %.4f %.4f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n" % 
				     	#		(i, j, thmin, thmax, self.f_bal[i,j], self.mean[i,j], 
			          	#		self.std_dev[i,j], self.ks_p_value[i,j], self.mu[i,j], self.sigma[i,j],
			          	#		self.chi2[i,j], self.mean_qsos[i,j]) )
						
						print 

		output_file.close()
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
			self.mu[i,j] = float(data[8])
			self.sigma[i,j] = float(data[9])
			self.chi2[i,j] = float(data[10])

		f.close()

		return 0



#def function_to_minimise(A, mu, sigma, )
def write_sim_to_file(sim, fname="simulation.out"):

	f = open(fname, "w")

	f.write("# SIMULATION OUTPUT\n")
	f.write("# i, j, thmin, thmax, f_bal, mean, std_dev, ks_p_value, delta_mu, chi2\n")

	# create an array to write
	for i, thmin in enumerate(sim.thetamin):
		for j,thmax in enumerate(sim.thetamax):
			f.write("%i %i %.4f %.4f %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n" % 
				     (i, j, thmin, thmax, sim.f_bal[i,j], sim.mean[i,j], 
			          sim.std_dev[i,j], sim.ks_p_value[i,j], sim.mean - sim.mean_qsos, sim.chi2) )

	f.close()



if __name__ == "__main__":

	# define angles to run sim over
	thetamin = np.arange(85,0,-5)
	thetamax = np.arange(90,0,-5)

	#thetamin = [84]
	#thetamax = [90]
	#thetamin = np.array([70])
	#thetamax = np.array([90])

	# get the data
	data = sub.get_sample("data/catalog.dat")

	# selection criteria
	# we need to be at a certain redshift and have O3 emission
	# then we select for BAL or no-BAL
	
	set_pretty()

	mode = sys.argv[1]

	source = sys.argv[2]

	LINE = "ew_o3"

	#SAVES = [(90,90), (20,40), (50,75), (30,80)]
	SAVES=None

	if source == "hst":

		# read the SHT COS catalog
		d_hst = sub.get_hst("data/hst_catalog.dat")

		# find sources which match between hst and SDSS
		matches, hst_matches, map_for_sdss, map_for_hst = sub.get_sdss_hst_matches(data, d_hst)
		
		# names of HST sources
		hst_names = np.loadtxt("data/hst_catalog.dat", usecols=(1,), unpack=True, dtype="string")

		# names of HST sources which match with SDSS
		hst_names_matches = hst_names[hst_matches]

		# names of BAL HST sources
		bal_names = np.loadtxt("HST/BALs", dtype="string")
		bal_flags = np.zeros(len(hst_names_matches))

		# find where the BALs are in HST
		for i, f in enumerate(hst_names_matches):
			select_bal_names = (f == bal_names)
			if np.sum(select_bal_names == 1): 
				bal_flags[i] = 1

		# turn the MAP into a flag of where BALs are in SDSS
		bal_use_flags = map_for_hst[bal_flags.astype(bool)]

		# pass this information to the selection class
		select = selection(data, hst_map = bal_use_flags)

	else:
		select = selection(data)


	if mode == "max" or mode == "faceon":

		print "here"

		# set up sim
		sim = simulation(thetamin, thetamax, data, select, line_string = LINE, mode=mode, source=source, saves = SAVES)

		sim.distribution = np.random.normal
		# run sim
		sim.run(select)
		#sim.run_gauss(select)

		write_sim_to_file(sim, fname="simulation_%s_%s_%s.out" % (LINE, mode, source))

		# make the contour plots
		#make_plots.plot_contour(sim)
		#make_plots.fourbytwo(sim)
		#make_plots.plot_contour2(sim)
		make_plots.p_max(sim)

	if mode == "nomax":

		# set up sim
		sim = simulation(thetamin, thetamax, data, select, line_string = LINE, mode=mode, source=source, saves = SAVES)

		# run sim
		sim.run(select)
		#sim.run_gauss(select)

		# make the contour plots
		#make_plots.plot_contour(sim)
		#write_sim_to_file(sim, fname="simulation_%s_%s_%s.out" % (LINE, mode, source))
		#make_plots.fourbytwo(sim)
		#make_plots.plot_contour2(sim)
		if mode == "faceon":
			make_plots.cont_faceon(sim)

	#elif mode == "read":

	elif mode == "hist":

		sel = selection(data)

		make_plots.make_hist2(data, sel)
		#make_plots.radio_hist(data, sel)
		#make_plots.radio_scatter(data, sel)
		#mixedCase
		#make_plots



