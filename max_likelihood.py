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




def get_dist(z):

	'''distance in pc from redshift'''

	return (z * 3e5 / 70.0) * 1e6


def emissivity_function(costheta):

	return costheta


def is_source_detected(costheta):

	dv = np.random.random()

	r = np.power(dv, 1.0/3.0)

	d = get_dist(r) * PARSEC

	flux = (1e46 * emissivity_function(costheta)) / (4.0*np.pi*d*d)

	d_detect = get_dist(1.0) * PARSEC
	detection_limit = 1e45 / (4.0*np.pi*d*d)

	#return (flux > detection_limit)
	return True


def get_mock_angles(THRESHOLD, NPTS, max_angle=None):

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


		costhetas[j] = costheta

		if (theta > THRESHOLD) and (theta < max_angle):
			bal_flags[j] = "b"
		else:
			bal_flags[j] = "q"

	return costhetas, bal_flags


def set_subplot_ticks():

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











# get the data
data = sub.get_sample("/Users/jmatthews/Documents/J_ApJS_194_45.tar/catalog.dat")

#some plotting modes
logbins = True
NORM  = True
lims = (0,150)
binsize = 5
labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA"]
set_pretty()
big_tick_labels(18)
colors=get_colors()


# selection criteria
# we need to be at a certain redshift and have O3 emission
# then we select for BAL or no-BAL
select = selection(data)

# these are the "threshold angles" in degrees
THRESHOLDS = [30.0, 45.0, 60.0, 75.0, 80.0]

### make the basic data plot
figure(figsize=(17,14))
subplot(2,3,1)
long_ticks()

# decide if we need logbins
bins = np.arange(lims[0],lims[1],binsize)

ews = data["ew_o3"]
if logbins: 
	ews = np.log10(ews)
	bins = np.arange(-2,4,0.1)

hist(ews[select.general*select.nonbal],bins=bins, facecolor=colors[0], alpha=0.5, log=True, label="non-BALs", normed=NORM, stacked=True)
hist(ews[select.general*select.mgbal],bins=bins, facecolor=colors[1], alpha=0.5, log=True, label="BALs", normed=NORM, stacked=True)
ylabel("Normalised Counts", fontsize=24)
set_subplot_ticks()

# now iterate over the thresold angles
for iang in range(len(THRESHOLDS)):

	subplot(2,3,iang+2)

	mock_data = data["ew_o3"][select_general]

	NPTS =	len(mock_data)

	THRESHOLD = THRESHOLDS[iang]

	# get angles above the threshold generated uniformly on sky
	# and apply to EW measurements
	costhetas, mock_bal_flags = get_mock_angles(THRESHOLD, NPTS, max_angle=None)

	mock_data = mock_data / emissivity_function(costhetas)
	select_mock_bals = (mock_bal_flags == "b")
	select_mock_nonbals = (mock_bal_flags == "q")

	bal_frac = float(np.sum(select_mock_bals)) / float(NPTS)
	

	ews = data["ew_o3"]

	if logbins: 
		mock_data = np.log10(mock_data)
		ews = np.log10(ews)

	# plot the original data
	hist(ews[select.general*select.nonbal],bins=bins, facecolor=colors[0], alpha=0.5, log=True, label="non-BALs", normed=NORM, stacked=True)
	hist(ews[select.general*select.mgbal],bins=bins, facecolor=colors[1], alpha=0.5, log=True, label="BALs", normed=NORM, stacked=True)

	# plot the mock data
	hist(mock_data[select_mock_nonbals], bins=bins, facecolor=colors[3], alpha=0.5, label=r"non-BAL Mock Data", normed=NORM, stacked=True)
	hist(mock_data[select_mock_bals], bins=bins, facecolor=colors[4], alpha=0.5, label=r"BAL Mock Data", normed=NORM, stacked=True)

	# decorate the figures a bit
	set_subplot_ticks()

	if iang == 2: ylabel("Normalised Counts", fontsize=24)
	
	ylimits = gca().get_ylim()
	xlimits = gca().get_xlim()
	text(0.6*xlimits[1], 0.4*ylimits[1],str(THRESHOLD), fontsize=24)

	float_legend(loc=2)

	'''
	If the K-S statistic is small or the p-value is high, 
	then we cannot reject the hypothesis that the distributions 
	of the two samples are the same.
	'''
	print "Angle ", THRESHOLD
	print '------------------'
	print stats.ks_2samp(mock_data[select_mock_bals], ews[select.general*select.nonbal])
	print stats.ks_2samp(mock_data[select_mock_nonbals], ews[select.general*select.mgbal])
	print bal_frac
	print ""




subplots_adjust(left=0.08,right=0.96)
savefig("mc_ew_exp_noselect.png")





