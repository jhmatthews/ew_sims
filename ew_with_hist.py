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
import seaborn as sns

import seaborn.apionly as sns

fig_size = (17,7)
#sns.set(context="paper", style="dark",rc={'text.usetex': True,'figure.figsize': fig_size, 'font.family': [u'serif'], 'font.serif':[u'Computer Modern'], 'axes.linewidth': 0.7, 'xtick.major.width': 0.7})
#sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

set_pretty()
big_tick_labels(18)
colors=get_colors()

# first get the data
columns = (13,14, 74, 87, 107)
bal_flag, r_flag, ews_o3, ews_mg2, ews_c4 = np.loadtxt("/Users/jmatthews/Documents/J_ApJS_194_45.tar/catalog.dat", unpack=True, usecols=columns)
select2 = (bal_flag == 0) 
select3 = (bal_flag >= 1) 

red_lims = [(0.35,0.83),(0.35,0.83),(1.45,2.28),(1.45,2.28)]

ews_to_do = [ews_o3, ews_mg2, ews_c4, ews_mg2]
n_to_do = len(ews_to_do)

lims = [(0,150),(0,200),(0,200)]
binsize = [5,5,5,5]
labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA"]

NORM  = True
THRESHOLDS = [30.0, 70.0]
logbins = True
# now make the histogram plot
figure(figsize=(20,7))

# loop over lines / subplots
for i in range(n_to_do):

	subplot(1,n_to_do, i+1)
	long_ticks()
	bins = np.arange(lims[i][0],lims[i][1],binsize[i])
	
	if logbins: bins = np.arange(-2,4,0.1)


	if i == 0: ylabel("Normalised Counts", fontsize=20)
	
	#ews = np.log10(ews_to_do[i])


	if logbins:
		ews = np.log10(ews_to_do[i])
	else:
		ews = ews_to_do[i]

	select1 = (ews>-3)

	select1 = (ews>0) 
	select2 = (bal_flag == 0) 
	select3 = (bal_flag >= 1) 

	hist(ews[select1*select2],bins=bins, facecolor=colors[0], alpha=0.7, log=True, label="non-BALs", normed=NORM, stacked=True)
	hist(ews[select1*select3],bins=bins, facecolor=colors[1], alpha=0.4, log=True, label="BALs", normed=NORM, stacked=True)

	ylimits = gca().get_ylim()
	text(0.4*lims[i][1], 0.6*ylimits[1],labels[i], fontsize=24)
	title(labels[i], fontsize=24)
	#ylim(0,0.06)
	xlabel(r"$\log [W_{\lambda}$ (\AA)]", fontsize=20)
	#xlim(lims[i][0],lims[i][1])

	text(0.25,4,r"$\mu_{non-BAL} = %.2f$\AA" % np.mean(10.0**ews[select1*select2]), fontsize=20)
	text(0.25,2,r"$\mu_{BAL} = %.2f$\AA" % np.mean(10.0**ews[select1*select3]), fontsize=20)

	if i == 0: 
		float_legend()

	# print stats.ks_2samp(ews[select1*select2], ews[select1*select3])

	# n1 = len(ews[select1*select2])
	# n2 = len(ews[select1*select3])
	# test = stats.ks_2samp(ews[select1*select2], ews[select1*select3])

	#print n1, n2
	#else:
	#gca().set_yticklabels([])

	data = ews[select1*select2]

	NPTS =	len(data)

	#mock_data = np.random.normal(mu, sigma, NPTS)

	for iang in range(len(THRESHOLDS)):

		if logbins:
			mock_data  = 10.0**data
		else:
			mock_data = data

		THRESHOLD = THRESHOLDS[iang]
		costhetas = np.zeros(NPTS)

		for j in range(NPTS):
			theta = 0
			while theta < THRESHOLD:

				costheta = np.random.random()
				theta = (np.arccos(costheta) * 180.0 / np.pi)

			costhetas[j] = costheta


		mock_data /= costhetas

		if logbins: mock_data = np.log10(mock_data)

		#hist(mock_data, bins=bins, facecolor=colors[2+iang], alpha=0.5, label=str(THRESHOLD), normed=NORM, stacked=True)
	#if logbins:
	#	(mu, sigma) = stats.norm.fit(ews[select1*select2])
	#	y = mlab.normpdf(bins, mu, sigma)
	#	l = plt.plot(bins, y, 'r--', linewidth=2)

	float_legend()
	xlim(0,2.1)
	#semilogy()

subplots_adjust(wspace = 0.15, left= 0.06, right=0.98, top=0.90)
savefig("ew_hist_qsos.png", dpi=300)

