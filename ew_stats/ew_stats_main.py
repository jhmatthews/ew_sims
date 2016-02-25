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
			#detected = True

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

	mu, sigma = p

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



def function_to_minimise(params, ew_o_quasars, costhetas):

	mu = params[0]
	sigma = params[1]
	sigma = 5.0

	ewstar = np.random.normal(loc=mu, scale=sigma, size=len(costhetas))

	ew_for_test = ewstar / costhetas

	'''
	If the K-S statistic is small or the p-value is high, 
	then we cannot reject the hypothesis that the distributions 
	of the two samples are the same.
	'''

	'''IMPROVE: is minimising 1/pval Ok?'''

	f_ew_for_test = histogram(ew_for_test, bins=np.arange(0,200,1))[0]
	f_ewo = histogram(ew_o_quasars, bins=np.arange(0,200,1))[0]

	select = (f_ewo > 5) * (f_ew_for_test > 5)

	df = np.sqrt(f_ewo)

	#chi2 = stats.chisquare(f_ewo, f_ew_for_test)

	chi2 = np.sum( (f_ewo[select] - f_ew_for_test[select])**2 / df[select] )

	return chi2



def check(ews, costhetas):

	mu = 5.0 + (20.0 * np.random.random(size=10000))
	sigma = (20.0 * np.random.random(size=10000))

	chi2_min = 1e100
	params_min = None

	for i, m in enumerate(mu):
		for j, s in enumerate(sigma):
			chi2 = function_to_minimise([m,s], ews, costhetas)

			if chi2 < chi2_min:
				chi2_min = chi2
				params_min = (m,s)
				print i, chi2, m, s

	return chi2_min, params_min



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
		redshift_lims_b = (3800.0 / 1550.0 - 1.0, 9200.0 / 2800.0 - 1.0)

		self.nonbal = (data["bal_flag"] == 0) 
		self.mgbal = (data["bal_flag"] == 2)
		self.bal = (data["bal_flag"] > 0)
		self.z = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])
		self.has_o3 = (data["ew_o3"] > 0)
		self.mass = None
		self.general = self.has_o3 * self.z
		self.a = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])
		self.b = (data["z"] > redshift_lims_b[0]) * (data["z"] < redshift_lims_b[1])



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
	def write_sim_to_file(self, fname="simulation.out"):

		f = open(fname, "w")

		f.write("# SIMULATION OUTPUT\n")
		f.write("# i, j, thmin, thmax, f_bal, mean, std_dev, ks_p_value\n")

		# create an array to write
		for i, thmin in enumerate(self.thetamin):
			for j,thmax in enumerate(self.thetamax):
				f.write("%i %i %.4f %.4f %8.4e %8.4e %8.4e %8.4e\n" % 
					     (i, j, thmin, thmax, self.f_bal[i,j], self.mean[i,j], self.std_dev[i,j], self.ks_p_value[i,j]) )

		f.close()


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


#def function_to_minimise(A, mu, sigma, )
def write_sim_to_file(sim, fname="simulation.out"):

	f = open(fname, "w")

	f.write("# SIMULATION OUTPUT\n")
	f.write("# i, j, thmin, thmax, f_bal, mean, std_dev, ks_p_value\n")

	# create an array to write
	for i, thmin in enumerate(sim.thetamin):
		for j,thmax in enumerate(sim.thetamax):
			f.write("%i %i %.4f %.4f %8.4e %8.4e %8.4e %8.4e\n" % 
				     (i, j, thmin, thmax, sim.f_bal[i,j], sim.mean[i,j], sim.std_dev[i,j], sim.ks_p_value[i,j]) )

	f.close()


def make_hist(data, select):

	strings = ["ew_o3", "ew_mg2", "ew_c4", "ew_mg2"]
	selects = [select.a, select.a, select.b, select.b]
	bal_selects = [select.mgbal, select.mgbal, select.bal, select.bal]
	logbins = True
	lims = [(0,150),(0,200),(0,200),(0,200)]
	labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA", r"Mg~\textsc{ii}~$2800$\AA"]
	NORM=True
	# now make the histogram plot
	set_pretty()
	figure(figsize=(20,7))

	for i in range(4):
		subplot(1,4, i+1)
		long_ticks()
		#bins = np.arange(lims[i][0],lims[i][1],binsize[i])

		if logbins: bins = np.arange(-2,4,0.1)


		if logbins:
			ews = np.log10(data[strings[i]])
		else:
			ews = ews_to_do[data[strings[i]]]
	

		hist(ews[selects[i]*select.nonbal],bins=bins, facecolor=colors[0], alpha=0.7, log=True, label="non-BALs", normed=NORM, stacked=True)
		hist(ews[selects[i]*bal_selects[i]],bins=bins, facecolor=colors[1], alpha=0.4, log=True, label="BALs", normed=NORM, stacked=True)

		if i == 0: ylabel("Normalised Counts", fontsize=20)

		xlim(0,3)
		ylim(1e-3,10)
		ylimits = gca().get_ylim()
		text(0.4*lims[i][1], 0.6*ylimits[1],labels[i], fontsize=20)
		title(labels[i], fontsize=24)
		#ylim(0,0.06)
		xlabel(r"$\log [W_{\lambda}$ (\AA)]", fontsize=20)
		#xlim(lims[i][0],lims[i][1])

		text(0.25,4,r"$\mu_{non-BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*select.nonbal]), fontsize=20)
		text(0.25,2,r"$\mu_{BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*bal_selects[i]]), fontsize=20)

		if i == 0:
			text(2,30,"Sample A, %i Mg BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=20)
		if i == 2:
			text(2,30,"Sample B, %i BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=20)

		if i == 0: 
			float_legend()


	subplots_adjust(wspace = 0.15, left= 0.06, right=0.98, top=0.85)
	savefig("ew_hist_qsos.png", dpi=300)



def make_2dhist(data, select):

	logbins = True
	lims = [(0,150),(0,200),(0,200),(0,200)]
	labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA", r"Mg~\textsc{ii}~$2800$\AA"]
	NORM=True
	# now make the histogram plot
	set_pretty()
	figure(figsize=(20,7))

	subplot(211)
	ews = np.log10(data["ew_c4"])
	hist(ews[selects[i]*select.nonbal],bins=bins, facecolor=colors[0], alpha=0.7, log=True, label="non-BALs", normed=NORM, stacked=True)
	hist(ews[selects[i]*bal_selects[i]],bins=bins, facecolor=colors[1], alpha=0.4, log=True, label="BALs", normed=NORM, stacked=True)





	for i in range(4):
		subplot(1,4, i+1)
		long_ticks()
		#bins = np.arange(lims[i][0],lims[i][1],binsize[i])

		if logbins: bins = np.arange(-2,4,0.1)


		if logbins:
			ews = np.log10(data[strings[i]])
		else:
			ews = ews_to_do[data[strings[i]]]
	

		hist(ews[selects[i]*select.nonbal],bins=bins, facecolor=colors[0], alpha=0.7, log=True, label="non-BALs", normed=NORM, stacked=True)
		hist(ews[selects[i]*bal_selects[i]],bins=bins, facecolor=colors[1], alpha=0.4, log=True, label="BALs", normed=NORM, stacked=True)

		if i == 0: ylabel("Normalised Counts", fontsize=20)

		xlim(0,3)
		ylim(1e-3,10)
		ylimits = gca().get_ylim()
		text(0.4*lims[i][1], 0.6*ylimits[1],labels[i], fontsize=20)
		title(labels[i], fontsize=24)
		#ylim(0,0.06)
		xlabel(r"$\log [W_{\lambda}$ (\AA)]", fontsize=20)
		#xlim(lims[i][0],lims[i][1])

		text(0.25,4,r"$\mu_{non-BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*select.nonbal]), fontsize=20)
		text(0.25,2,r"$\mu_{BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*bal_selects[i]]), fontsize=20)

		if i == 0:
			text(2,30,"Sample A, %i Mg BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=20)
		if i == 2:
			text(2,30,"Sample B, %i BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=20)

		if i == 0: 
			float_legend()


	subplots_adjust(wspace = 0.15, left= 0.06, right=0.98, top=0.85)
	savefig("ew_hist_qsos.png", dpi=300)



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

	mode = sys.argv[1]

	if mode == "sim":

		# set up sim
		sim = simulation(thetamin, thetamax, data)

		# run sim
		sim.run(select)
		#sim.run_gauss(select)

		# make the contour plots
		make_plots.plot_contour(sim)
		make_plots.plot_contour2(sim)

	#elif mode == "read":

	elif mode == "hist":

		make_hist(data, select)


