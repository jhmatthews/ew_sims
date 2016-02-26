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

def individual_histogram(line, thmin, thmax, ew_o, ew_mock, chi2, bins):

	colors = get_colors()

	fig1=figure()

	frame1=fig1.add_axes((.1,.3,.8,.6))
	NORM = False
	n1, bins1, patches1 = hist(ew_o, normed=NORM, facecolor=colors[0], alpha=0.5, bins=bins, label="EW$_O$")
	n2, bins2, patches2 = hist(ew_mock, normed=NORM, facecolor=colors[1], alpha=0.5, bins=bins, label="EW$_{fit}$")
	float_legend()
	ylabel("$N$", fontsize=20)
	text(110,800,r"$\chi^2/dof=%.2f$" % chi2, fontsize=20)
	chi = (n1 - n2) / np.sqrt(n1)
	frame1.set_xticklabels([])
	title(r"$\theta_{min} = %i, \theta_{max} = %i$" % (thmin, thmax), fontsize=20)
	
	frame2=fig1.add_axes((.1,.1,.8,.2))
	plot(0.5*(bins1[1:] + bins1[:-1]), chi, linewidth=2, c="k")
	ylabel("$\chi$", fontsize=20)
	xlabel("EW (\AA)", fontsize=20)
	savefig("hists_%s/hist_%i_%i_%s.png" % (line,thmin, thmax, line), dpi=300)
	clf()

	return 0


def bal_histogram(line, thmin, thmax, ew_o, ew_mock, pval, f_bal, bins):

	colors = get_colors()

	fig1=figure()

	NORM = True
	n1, bins1, patches1 = hist(ew_o, normed=NORM, facecolor=colors[0], alpha=0.5, bins=bins, label="EW$_O$")
	n2, bins2, patches2 = hist(ew_mock, normed=NORM, facecolor=colors[1], alpha=0.5, bins=bins, label="EW$_{fit}$")
	float_legend()
	ylabel("$N/N_{tot}$", fontsize=20)
	ylim = gca().get_ylim()
	text(80,0.6*np.sum(ylim),r"$p_{KS}=%8.4e$" % (pval), fontsize=20)
	text(80,0.45*np.sum(ylim),r"$f_{BAL}=%.2f$" % (f_bal), fontsize=20)
	


	title(r"BALs, $\theta_{min} = %i, \theta_{max} = %i$" % (thmin, thmax), fontsize=20)
	xlabel("EW (\AA)", fontsize=20)

	savefig("hists_%s/balhist_%i_%i_%s.png" % (line,thmin, thmax, line), dpi=300)
	clf()

	return 0




def plot_contour(sim):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(20,7))
	subplot(1,3,1)

	sim.ks_p_value = np.ma.masked_array(sim.ks_p_value, mask=(sim.ks_p_value == 0) )
	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )
	sim.chi2 = np.ma.masked_array(sim.chi2, mask=(sim.chi2 == 0) )

	contourf(sim.thetamax, sim.thetamin, np.log10(sim.ks_p_value), extend='both', levels=np.arange(-5,0,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, np.log10(ks_test), levels=np.log10(np.array([0.001,0.05,0.32])), c="k")
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title('''$\log~(p$-value$)$. If the p-value is high, 
	then we cannot reject the hypothesis that the distributions 
	of the two samples are the same''', fontsize=14)

	subplot(1,3,2)
	contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title("$f_{BAL}$", fontsize=20)

	subplot(1,3,3)
	contourf(sim.thetamax, sim.thetamin, sim.chi2, 
		     extend='both', levels=np.arange(0,50,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\chi^2 / dof$", fontsize=20)

	subplots_adjust(left=0.1,right=0.97,top=0.80)
	savefig("contour_%s.png" % sim.line_string, dpi=200)

	return 0


def plot_contour2(sim):

	'''
	sim 
	instance of class simulation 
	'''

	d_use = sim.data[sim.line_string][sim.s_bal]

	mu = np.mean(d_use)
	sigma = np.std(d_use)

	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.std_dev = np.ma.masked_array(sim.std_dev, mask=(sim.std_dev == 0) )

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(13,7))
	subplot(1,2,1)

	contourf(sim.thetamax, sim.thetamin, sim.mean - mu, extend="max", cmap=cm_use, levels=np.arange(-10,50,1))
	colorbar()

	CS4 =contour(sim.thetamax, sim.thetamin, sim.mean - mu, colors=('w',), linewidths=(2,),levels=np.arange(-10,60,10))
	clabel(CS4, fmt='%2i', colors='w', fontsize=14)

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\Delta \mu = \mu_{EW,mock} - \mu_{EW,BALs}$", fontsize=20)

	subplot(1,2,2)
	contourf(sim.thetamax, sim.thetamin, cp - sigma, extend="max", cmap=cm_use, levels=np.arange(-10,50,1))
	colorbar()

	CS4 = contour(sim.thetamax, sim.thetamin, sim.std_dev - sigma,  colors=('w',), linewidths=(2,), levels=np.arange(-10,60,10))
	#contour(thetamax, thetamins, f_bal)
	clabel(CS4, fmt='%2i', colors='w', fontsize=14)
	
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\Delta \sigma = \sigma_{EW,mock} - \sigma_{EW,BALs}$", fontsize=20)

	subplots_adjust(left=0.1,right=0.97,top=0.93)
	savefig("contour2_%s.png" % sim.line_string, dpi=200)

	return 0




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

def set_subplot_ticks():
	'''
	plotting hack
	'''

	xlabel(r"$W_{\lambda}$ (\AA)", fontsize=24)
	xlim(-1,4)
	ylim(1e-3,2)
	gca().set_xticklabels(["0.1","1","10","$10^2$","$10^3$","$10^4$"])

	return
