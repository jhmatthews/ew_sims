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

def individual_histogram(sim, thmin, thmax, ew_o, ew_mock, chi2, bins):

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
	plot(0.5*(bins1[1:] + bins1[:-1]), chi, linewidth=2, c="k", ls="steps")
	ylabel("$\chi$", fontsize=20)
	xlabel("EW (\AA)", fontsize=20)
	#savefig("hists_%s/hist_%i_%i_%s.png" % (line,thmin, thmax, line), dpi=300)
	#savefig("examples_%s_%s_%s/hist_%i_%i_%s.png" % (sim.line_string, sim.mode, sim.source, thmin, thmax, sim.line_string), dpi=300)
	savefig("hist_quasar_fit.png", dpi=300)
	clf()

	return 0

def f_fit(x, alpha, norm):

	'''
	power law of form norm norm * (x**alpha)
	used for fitting continuum 
	'''

	return norm * (x**alpha)

import scipy.optimize as opt


def individual_scatter(sim, thmin, thmax, ew_o, ew_mock, chi2, bins):

	colors = get_colors()

	fig1=figure()

	frame1=fig1.add_axes((.1,.3,.8,.6))
	NORM = False
	n1, bins1 = histogram(ew_o, normed=NORM, bins=bins)
	nerr = np.sqrt(n1)
	bin_centres = 0.5* (bins1[1:] + bins1[:-1])

	errorbar(bin_centres, n1, yerr=nerr, fmt=None, label="EW$_O$", capsize=0, elinewidth=1.5, capthick=1.5, ecolor=colors[2])


	n2, bins2 = histogram(ew_mock, normed=NORM, bins=bins)
	plot(bins1[:-1], n2, c=colors[0], linewidth=2, ls="steps")

	coeffs = opt.curve_fit(f_fit, bin_centres[30:], n1[30:])

	eww = np.arange(1,150,0.1)
	plot(eww, f_fit(eww, coeffs[0][0], coeffs[0][1]), linestyle="--", linewidth=2, c="k")
	print coeffs[0]


	loglog()
	xlim(2,150)
	ylim(0.5,1e4)
	
	float_legend()
	ylabel("$N$", fontsize=20)
	text(10,10,r"$\chi^2/dof=%.2f$" % chi2, fontsize=20)
	chi = (n1 - n2) / np.sqrt(n1)
	frame1.set_xticklabels([])
	#title(r"$\theta_{min} = %i, \theta_{max} = %i$" % (thmin, thmax), fontsize=20)
	
	frame2=fig1.add_axes((.1,.1,.8,.2))

	plot(0.5*(bins1[1:] + bins1[:-1]), chi, linewidth=2, c="k", ls="steps")
	hlines([0.0], 0,160, linestyle="--")
	ylabel("$\chi$", fontsize=20)
	xlabel("EW (\AA)", fontsize=20)
	semilogx()
	xlim(2,150)
	#savefig("hists_%s/hist_%i_%i_%s.png" % (line,thmin, thmax, line), dpi=300)
	#savefig("examples_%s_%s_%s/hist_%i_%i_%s.png" % (sim.line_string, sim.mode, sim.source, thmin, thmax, sim.line_string), dpi=300)
	savefig("quasar_fit.png", dpi=300)
	clf()

	return 0


def bal_histogram(sim, thmin, thmax, ew_o, ew_mock, mu, f_bal, bins):

	colors = get_colors()

	fig1=figure()

	NORM = True
	n1, bins1, patches1 = hist(ew_o, normed=NORM, facecolor=colors[0], alpha=0.5, bins=bins, label="EW$_O$")
	n2, bins2, patches2 = hist(ew_mock, normed=NORM, facecolor=colors[1], alpha=0.5, bins=bins, label="EW$_{fit}$")
	float_legend()
	ylabel("$N/N_{tot}$", fontsize=20)
	ylim = gca().get_ylim()
	text(80,0.6*np.sum(ylim),r"$\mu=%8.4e$" % (mu), fontsize=20)
	#text(80,0.6*np.sum(ylim),r"$p_{KS}=%8.4e$" % (pval), fontsize=20)
	text(80,0.45*np.sum(ylim),r"$f_{BAL}=%.2f$" % (f_bal), fontsize=20)
	


	title(r"BALs, $\theta_{min} = %i, \theta_{max} = %i$" % (thmin, thmax), fontsize=20)
	xlabel("EW (\AA)", fontsize=20)

	#savefig("hists_%s/balhist_%i_%i_%s.png" % (line,thmin, thmax, line), dpi=300)
	savefig("examples_%s_%s_%s/balhist_%i_%i_%s.png" % (sim.line_string, sim.mode, sim.source, thmin, thmax, sim.line_string), dpi=300)
	clf()

	return 0




def twobytwo_contour(sim):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(10,7))
	subplot(2,2,1)

	sim.ks_p_value = np.ma.masked_array(sim.ks_p_value, mask=(sim.ks_p_value == 0) )
	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )
	sim.chi2 = np.ma.masked_array(sim.chi2, mask=(sim.chi2 == 0) )

	contourf(sim.thetamax, sim.thetamin, np.log10(sim.ks_p_value), extend='both', levels=np.arange(-3,0,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, np.log10(ks_test), levels=np.log10(np.array([0.001,0.05,0.32])), c="k")
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	#title('''$\log~(p$-value$)$. If the p-value is high, 
	#then we cannot reject the hypothesis that the distributions 
	#of the two samples are the same''', fontsize=14)

	subplot(2,2,2)
	contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,0.6,0.02), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title("$f_{BAL}$", fontsize=20)

	subplot(2,2,3)
	contourf(sim.thetamax, sim.thetamin, sim.chi2, 
		     extend='both', levels=np.arange(0,10,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\chi^2 / dof$", fontsize=20)


	subplot(2,2,4)
	d_use = sim.data[sim.line_string][sim.s_bals]

	mu = np.mean(d_use)
	sigma = np.std(d_use)

	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.std_dev = np.ma.masked_array(sim.std_dev, mask=(sim.std_dev == 0) )

	contourf(sim.thetamax, sim.thetamin, sim.mean - mu, extend="max", cmap=cm_use, levels=np.arange(-50,50,1))
	colorbar()

	CS4 =contour(sim.thetamax, sim.thetamin, sim.mean - mu, colors=('k',), linewidths=(2,),levels=np.arange(-50,60,10))
	clabel(CD4, fmt='%2i', colors='k', fontsize=14)

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\Delta \mu = \mu_{EW,mock} - \mu_{EW,BALs}$", fontsize=20)


	#subplots_adjust(left=0.1,right=0.97,top=0.80)
	savefig("twobytwo_%s_%s.png" % (sim.line_string, sim.mode), dpi=200)

	return 0





def plot_contour(sim):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(10,7))
	subplot(1,3,1)

	sim.ks_p_value = np.ma.masked_array(sim.ks_p_value, mask=(sim.ks_p_value == 0) )
	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )
	sim.chi2 = np.ma.masked_array(sim.chi2, mask=(sim.chi2 == 0) )

	contourf(sim.thetamax, sim.thetamin, np.log10(sim.ks_p_value), extend='both', levels=np.arange(-5,0,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, np.log10(ks_test), levels=np.log10(np.array([0.001,0.05,0.32])), c="k")
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	#title('''$\log~(p$-value$)$. If the p-value is high, 
	#then we cannot reject the hypothesis that the distributions 
	#of the two samples are the same''', fontsize=14)

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
	savefig("contour_%s_%s.png" % (sim.line_string, sim.mode), dpi=200)

	return 0


def plot_contour2(sim):

	'''
	sim 
	instance of class simulation 
	'''

	d_use = sim.data[sim.line_string][sim.s_bals]

	mu = np.mean(d_use)
	sigma = np.std(d_use)

	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.std_dev = np.ma.masked_array(sim.std_dev, mask=(sim.std_dev == 0) )

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(13,7))
	subplot(1,2,1)

	contourf(sim.thetamax, sim.thetamin, sim.mean - mu, extend="max", cmap=cm_use, levels=np.arange(-50,50,1))
	colorbar()

	CS4 =contour(sim.thetamax, sim.thetamin, sim.mean - mu, colors=('w',), linewidths=(2,),levels=np.arange(-50,60,10))
	clabel(CS4, fmt='%2i', colors='w', fontsize=14)

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\Delta \mu = \mu_{EW,mock} - \mu_{EW,BALs}$", fontsize=20)

	subplot(1,2,2)
	contourf(sim.thetamax, sim.thetamin, sim.std_dev - sigma, extend="max", cmap=cm_use, levels=np.arange(-50,50,1))
	colorbar()

	CS4 = contour(sim.thetamax, sim.thetamin, sim.std_dev - sigma,  colors=('w',), linewidths=(2,), levels=np.arange(-50,60,10))
	#contour(thetamax, thetamins, f_bal)
	clabel(CS4, fmt='%2i', colors='w', fontsize=14)
	
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title(r"$\Delta \sigma = \sigma_{EW,mock} - \sigma_{EW,BALs}$", fontsize=20)

	subplots_adjust(left=0.1,right=0.97,top=0.93)
	savefig("contour2_%s_%s.png" % (sim.line_string, sim.mode), dpi=200)

	return 0


def make_hist(data, select):

	selects = [select.a, select.a, select.b, select.b, select.a]
	bal_selects = [select.mgbal, select.mgbal, select.bal, select.bal, select.bal]
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
	

		hist(ews[selects[i]*select.nonbal],bins=bins, ec=colors[0], lw = 3, facecolor="None", alpha=1, log=True, label="non-BALs", normed=NORM, stacked=True)
		hist(ews[selects[i]*bal_selects[i]],bins=bins, ec=colors[2], lw = 3, facecolor="None", alpha=1, log=True, label="BALs", normed=NORM, stacked=True)

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


def radio_hist(data, select):

	strings = ["fwhm_hb", "ew_mg2"]
	selects = [select.a*select.lobe, select.a*select.core]
	bal_selects = [select.mgbal, select.mgbal, select.bal, select.bal]
	logbins = False
	lims = [(0,150),(0,200),(0,200),(0,200)]
	labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA", r"Mg~\textsc{ii}~$2800$\AA"]
	NORM=True
	# now make the histogram plot
	set_pretty()
	colors=[get_colors()[0], get_colors()[4]]
	figure()
	big_tick_labels(16)


	#subplot(2,2, i+1)
	long_ticks()
	#bins = np.arange(lims[i][0],lims[i][1],binsize[i])

	if logbins: 
		bins = np.arange(-2,4,0.1)
	else:
		bins = np.arange(0,15000,500)


	if logbins:
		ews = np.log10(data[strings[0]])
	else:
		ews = data[strings[0]]
	

	hist(ews[selects[0]*select.nonbal],bins=bins, facecolor=colors[0], alpha=0.7, log=True, label="Lobe", normed=NORM, stacked=True)
	hist(ews[selects[1]*select.nonbal],bins=bins, facecolor=colors[1], alpha=0.4, log=True, label="Core", normed=NORM, stacked=True)

		#if i == 0 or i == 2: 
	ylabel("Normalised Counts", fontsize=20)
		#else:
		#	gca().set_yticklabels([])

	#xlim(0,2.9)
	#ylim(1e-3,10)
	ylimits = gca().get_ylim()
	#text(0.4*lims[0][1], 0.6*ylimits[1],labels[0], fontsize=20)
	#title(labels[0], fontsize=20)
	#ylim(0,0.06)
	xlabel(r"$\log$ [EW (\AA)]", fontsize=20)
	#xlim(lims[i][0],lims[i][1])

	#text(0.25,4,r"$\mu_{non-BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*select.nonbal]), fontsize=20)
	#text(0.25,2,r"$\mu_{BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*bal_selects[i]]), fontsize=20)
	#if i == 0: 
	float_legend()

	long_ticks()
	#subplots_adjust(hspace = 0.35, wspace=0.05, left= 0.1, right=0.98, top=0.92, bottom=0.06)
	savefig("ew_radio.png", dpi=300)




def radio_scatter(data, select):

	strings = ["fwhm_hb", "ew_o3"]
	selects = [select.a*select.lobe*(data[strings[0]] > 0), select.a*select.core*(data[strings[0]] > 0)]
	bal_selects = [select.mgbal, select.mgbal, select.bal, select.bal]
	# now make the histogram plot
	set_pretty()
	cols=get_colors()
	figure()
	big_tick_labels(16)


	#subplot(2,2, i+1)
	long_ticks()
	#bins = np.arange(lims[i][0],lims[i][1],binsize[i])



	ews = np.log10(data[strings[1]])
	fws = np.log10(data[strings[0]])

	scatter(fws[selects[0]*select.nonbal], ews[selects[0]*select.nonbal],
		    edgecolors="k", facecolors="None", label="Lobe",alpha=0.5)
	scatter(fws[selects[1]*select.nonbal], ews[selects[1]*select.nonbal],
		    edgecolors=cols[1], facecolors="None", label="Core",alpha=0.5)
	

	scatter(fws[selects[0]*select.mgbal], ews[selects[0]*select.mgbal],
		    c=cols[2], marker="s", label="Lobe BAL", s=80)
	#scatter(fws[selects[1]*select.mgbal], ews[selects[1]*select.mgbal],
	#	    c=colors[4], marker="s", label="Core BAL", s=80)

		#if i == 0 or i == 2: 
	#ylabel("Normalised Counts", fontsize=20)
		#else:
		#	gca().set_yticklabels([])

	#xlim(0,2.9)
	#ylim(1e-3,10)
	#ylimits = gca().get_ylim()
	#text(0.4*lims[0][1], 0.6*ylimits[1],labels[0], fontsize=20)
	#title(labels[0], fontsize=20)
	#ylim(0,0.06)
	#xlabel(r"$\log$ [EW (\AA)]", fontsize=20)
	#xlim(lims[i][0],lims[i][1])

	#text(0.25,4,r"$\mu_{non-BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*select.nonbal]), fontsize=20)
	#text(0.25,2,r"$\mu_{BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*bal_selects[i]]), fontsize=20)
	#if i == 0: 
	pretty_legend()
	xlabel(r"EW[O~\textsc{iii} (\AA)]", fontsize=20)
	ylabel(r"FWHM[H$\beta$] (km~s$^{-1}$)", fontsize=20)

	long_ticks()
	#subplots_adjust(hspace = 0.35, wspace=0.05, left= 0.1, right=0.98, top=0.92, bottom=0.06)
	savefig("scatter_radio.png", dpi=300)





def make_hist2(data, select):

	strings = ["ew_o3", "ew_mg2", "ew_c4", "ew_mg2"]
	selects = [select.a, select.a, select.b, select.b]
	bal_selects = [select.mgbal, select.mgbal, select.bal, select.bal]
	logbins = True
	lims = [(0,150),(0,200),(0,200),(0,200)]
	labels=[r"[O~\textsc{iii}]~$5007$\AA", r"Mg~\textsc{ii}~$2800$\AA", r"C~\textsc{iv}~$1550$\AA", r"Mg~\textsc{ii}~$2800$\AA"]
	NORM=True
	# now make the histogram plot
	set_pretty()
	colors=get_colors()
	figure(figsize=(10,14))
	big_tick_labels(16)

	for i in range(4):
		subplot(2,2, i+1)
		long_ticks()
		#bins = np.arange(lims[i][0],lims[i][1],binsize[i])

		if logbins: bins = np.arange(-2,4,0.1)


		if logbins:
			ews = np.log10(data[strings[i]])
		else:
			ews = ews_to_do[data[strings[i]]]
	

		hist(ews[selects[i]*select.nonbal],bins=bins, facecolor=colors[0], lw = 3, ec="None", alpha=0.9, log=True, label="non-BALs", normed=NORM, stacked=True)
		hist(ews[selects[i]*bal_selects[i]],bins=bins, facecolor=colors[2], lw = 3, ec=colors[2], alpha=0.5, log=True, label="BALs", normed=NORM, stacked=True)

		if i == 0 or i == 2: 
			ylabel("Normalised Counts", fontsize=20)
		else:
			gca().set_yticklabels([])

		xlim(0,2.9)
		ylim(1e-3,10)
		ylimits = gca().get_ylim()
		text(0.4*lims[i][1], 0.6*ylimits[1],labels[i], fontsize=20)
		title(labels[i], fontsize=20)
		#ylim(0,0.06)
		xlabel(r"$\log$ [EW (\AA)]", fontsize=20)
		#xlim(lims[i][0],lims[i][1])

		text(0.25,4,r"$\mu_{non-BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*select.nonbal]), fontsize=20)
		text(0.25,2,r"$\mu_{BAL} = %.2f$\AA" % np.mean(10.0**ews[selects[i]*bal_selects[i]]), fontsize=20)

		if i == 0:
			text(1.6,30,"Sample A, %i Mg BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=24)
		if i == 2:
			text(1.6,30,"Sample B, %i BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=24)

		if i == 0: 
			float_legend()

	long_ticks()
	subplots_adjust(hspace = 0.35, wspace=0.05, left= 0.1, right=0.98, top=0.92, bottom=0.06)
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
			text(2,30,"Sample A, %i Mg BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=24)
		if i == 2:
			text(2,30,"Sample B, %i BALs, %i non-BALs." % ( np.sum(selects[i]*bal_selects[i]), np.sum(selects[i]*select.nonbal) ) , fontsize=24)

		if i == 0: 
			float_legend()


	subplots_adjust(wspace = 0.15, left= 0.06, right=0.98, top=0.80)
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

def fourbytwo(sim, mode="minmax"):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(10,12))
	subplot(4,2,1)

	sim.ks_p_value = np.ma.masked_array(sim.ks_p_value, mask=(sim.ks_p_value == 0) )
	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )
	sim.chi2 = np.ma.masked_array(sim.chi2, mask=(sim.chi2 == 0) )
	sim.mu = np.ma.masked_array(sim.mu, mask=(sim.mu == 0) )
	sim.sigma = np.ma.masked_array(sim.sigma, mask=(sim.sigma == 0) )
	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.std_dev = np.ma.masked_array(sim.std_dev, mask=(sim.std_dev == 0) )

	# if mode == "minmax":
	# 	x = thetamin 
	# 	y = thetamax

	contourf(sim.thetamax, sim.thetamin, np.log10(sim.ks_p_value), extend='both', levels=np.arange(-3,0,0.1), cmap=cm_use)
	colorbar()

	CS4 = contour(sim.thetamax, sim.thetamin, sim.ks_p_value, levels=[0.1,0.05])
	clabel(CS4, fmt='%.2f', colors='w', fontsize=14)
	#contour(thetamax, thetamins, np.log10(ks_test), levels=np.log10(np.array([0.001,0.05,0.32])), c="k")
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(30,70,r'$\log~(p_{KS}$-value$)$', fontsize=18)
	gca().set_xticklabels([])

	subplot(4,2,2)
	contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,0.6,0.02), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(30,70,"$f_{BAL}$", fontsize=20)
	gca().set_xticklabels([])

	subplot(4,2,3)
	contourf(sim.thetamax, sim.thetamin, sim.chi2, 
		     extend='both', levels=np.arange(0,30,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(30,70,r"$\chi^2 / dof$", fontsize=20)
	gca().set_xticklabels([])

	subplot(4,2,4)

	contourf(sim.thetamax, sim.thetamin, sim.mean, extend="max", cmap=cm_use, levels=np.arange(0,100,1))
	colorbar()

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(10,70,r"$\mu_{EW,mock}$", fontsize=20)
	gca().set_xticklabels([])


	subplot(4,2,5)
	d_use = sim.data[sim.line_string][sim.s_bals]

	mu = np.mean(d_use)
	sigma = np.std(d_use)

	contourf(sim.thetamax, sim.thetamin, sim.mean - mu, extend="max", cmap=cm_use, levels=np.arange(-50,50,1))
	colorbar()

	CS4 =contour(sim.thetamax, sim.thetamin, sim.mean - mu, colors=('w',), linewidths=(2,),levels=np.arange(-50,60,10))
	clabel(CS4, fmt='%2i', colors='w', fontsize=14)

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(10,70,r"$\mu_{EW,mock} - \mu_{EW,BALs}$", fontsize=15)
	gca().set_xticklabels([])

	subplot(4,2,6)
	d_use = sim.data[sim.line_string][sim.s_nonbals]

	mu = np.mean(d_use)
	sigma = np.std(d_use)

	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.std_dev = np.ma.masked_array(sim.std_dev, mask=(sim.std_dev == 0) )

	contourf(sim.thetamax, sim.thetamin, sim.mean - mu, extend="max", cmap=cm_use, levels=np.arange(-50,50,1))
	colorbar()

	CS4 =contour(sim.thetamax, sim.thetamin, sim.mean - mu, colors=('w',), linewidths=(2,),levels=np.arange(-50,60,10))
	clabel(CS4, fmt='%2i', colors='w', fontsize=14)

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(10,70,r"$\mu_{EW,mock} - \mu_{EW,non-BALs}$", fontsize=15)
	gca().set_xticklabels([])



	subplot(4,2,7)
	contourf(sim.thetamax, sim.thetamin, sim.mu, extend="max", cmap=cm_use, levels=np.arange(0,20,0.4))
	colorbar()

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(30,70,r"$\mu^*$", fontsize=20)

	subplot(4,2,8)
	contourf(sim.thetamax, sim.thetamin, sim.sigma, extend="max", cmap=cm_use, levels=np.arange(0,10,0.4))
	colorbar()

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(30,70,r"$\sigma^*$", fontsize=20)


	subplots_adjust(left=0.1,right=0.97,top=0.97, bottom=0.06, wspace=0.1, hspace=0.02)
	savefig("four_%s_%s_%s.png" % (sim.line_string, sim.mode, sim.source), dpi=200)

	return 0





def interpsmooth(x, y, val):

	index = np.indices(val.shape)

	#xx = x[index[0]]
	#yy = y[index[0]]

	interp_func = interp2d(xx, yy, val, kind='cubic')

	return interp_func





def p_max(sim, mode="minmax"):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(10,10))

	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )
	sim.ks_p_value = np.ma.masked_array(sim.ks_p_value, mask=(sim.ks_p_value == 0) )
	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.mean_qsos = np.ma.masked_array(sim.mean_qsos, mask=(sim.mean_qsos == 0) )
	sim.chi2 = np.ma.masked_array(sim.chi2, mask=(sim.chi2 == 0) )

	subplot(2,2,1)
	#contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
	pcolormesh(sim.thetamax, sim.thetamin, np.log10(sim.ks_p_value), cmap=cm_use, vmin=-4,vmax=0)
	colorbar()

	#contour(sim.thetamax, sim.thetamin, sim.f_bal, colors=('k',), linewidths=(2,), levels=[0.14,0.26,0.6])
	ylabel(r"$\theta_{min}$", fontsize=20)
	#xlabel(r"$\theta_{max}$", fontsize=20)
	#text(30,70,"$f_{BAL}$", fontsize=20)
	gca().set_xticklabels([])
	text(10,70,r"$\log(p_{KS})$", fontsize=20)
	ylim(5,90)


	subplot(2,2,2)
	#contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
	pcolormesh(sim.thetamax, sim.thetamin, sim.chi2, cmap=cm_use, vmin=0,vmax=50)
	colorbar()
	#contour(sim.thetamax, sim.thetamin, sim.f_bal, colors=('k',), linewidths=(2,), levels=[0.14,0.26,0.6])
	#ylabel(r"$\theta_{min}$", fontsize=20)
	#xlabel(r"$\theta_{max}$", fontsize=20)
	#text(30,70,"$f_{BAL}$", fontsize=20)
	gca().set_yticklabels([])
	gca().set_xticklabels([])
	text(10,70,r"$\chi^2 / dof$", fontsize=20)
	#gca().set_xticklabels([])
	ylim(5,90)

	subplot(2,2,3)
	#contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
	pcolormesh(sim.thetamax, sim.thetamin, sim.f_bal, cmap=cm_use, vmin=0,vmax=0.99)
	colorbar()
	#contour(sim.thetamax, sim.thetamin, sim.f_bal, colors=('k',), linewidths=(2,), levels=[0.14,0.26,0.6])
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	#text(30,70,"$f_{BAL}$", fontsize=20)
	text(10,70,r"$f_{BAL}$", fontsize=20)
	ylim(5,90)


	subplot(2,2,4)
	#contourf(sim.thetamax, sim.thetamin, sim.mean - sim.mean_qsos, extend="max", cmap=cm_use, levels=np.arange(0,100,1))
	pcolormesh(sim.thetamax, sim.thetamin, sim.mean - sim.mean_qsos, cmap=cm_use, vmin=0,vmax=29)
	colorbar()
	gca().set_yticklabels([])
	#ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(10,70,r"$\Delta \mu_{EW}$", fontsize=20)
	#gca().set_xticklabels([])
	ylim(5,90)



	subplots_adjust(left=0.1,right=0.97,top=0.97, bottom=0.06, wspace=0.05, hspace=0.05)
	savefig("mesh4_%s_%s_%s.png" % (sim.line_string, sim.mode, sim.source), dpi=200)

	return 0









def cont_faceon(sim, mode="minmax"):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(7,10))

	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )
	sim.mean = np.ma.masked_array(sim.mean, mask=(sim.mean == 0) )
	sim.mean_qsos = np.ma.masked_array(sim.mean_qsos, mask=(sim.mean_qsos == 0) )


	subplot(2,1,1)
	#contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
	pcolormesh(sim.thetamax, sim.thetamin, sim.f_bal, cmap=cm_use, vmin=0,vmax=1)
	colorbar()

	#contour(sim.thetamax, sim.thetamin, sim.f_bal, colors=('k',), linewidths=(2,), levels=[0.14,0.26,0.6])


	ylabel(r"$\theta_{min}$", fontsize=20)
	#xlabel(r"$\theta_{max}$", fontsize=20)
	#text(30,70,"$f_{BAL}$", fontsize=20)
	gca().set_xticklabels([])
	text(10,70,r"$f_{BAL}$", fontsize=20)
	ylim(5,90)


	subplot(2,1,2)

	#contourf(sim.thetamax, sim.thetamin, sim.mean - sim.mean_qsos, extend="max", cmap=cm_use, levels=np.arange(0,100,1))
	pcolormesh(sim.thetamax, sim.thetamin, sim.mean - sim.mean_qsos, cmap=cm_use, vmin=0,vmax=60)
	colorbar()

	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	text(10,70,r"$\Delta \mu_{EW}$", fontsize=20)
	#gca().set_xticklabels([])
	ylim(5,90)



	subplots_adjust(left=0.1,right=0.97,top=0.97, bottom=0.06, wspace=0.1, hspace=0.02)
	savefig("faceon_%s_%s_%s.png" % (sim.line_string, sim.mode, sim.source), dpi=200)

	return 0

