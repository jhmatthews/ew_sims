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


set_pretty()
rcParams["text.usetex"] = "True"
big_tick_labels(18)
colors=get_colors()
NORM  = True
THRESHOLDS = [30.0, 70.0]
logbins = False
LOG = True
# now make the histogram plot

data = sub.get_sample("/Users/jmatthews/Documents/J_ApJS_194_45.tar/catalog.dat")
to_plot = np.array(["Lbol", "L5100", "L3000", "z", "l_o3", "ew_o3", "ew_mg2", "mbh_hbeta"])
usebins = [np.arange(44,47,0.2),np.arange(44,47,0.2), np.arange(44,47,0.2),np.arange(0,7,0.1),
np.arange(40,44,0.2), np.arange(-1,3,0.1), np.arange(-1,3,0.1), np.arange(7.5,10,0.1)]
data_labels = sub.make_label_map()

logs = [False, False, False, False, False, True, True, False]

select_mass = False
select_z = False

redshift_lims = (3800.0 / 2800.0 - 1.0, 9200.0 / 5007.0 - 1.0)


def make_hist():
	figure(figsize=(12,18))
	for i in range(len(to_plot)):

		subplot(4,2,i+1)

		data_to_plot = data[to_plot[i]]

		select_null = (data_to_plot > 0) 

		if logs[i]:
			data_to_plot = np.log10(data_to_plot)

		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		if select_mass: select_null *= (data["mbh_hbeta"] > 8.5) * (data["mbh_hbeta"] < 9.5)
		select1 = (data["bal_flag"] == 0) 
		select2 = (data["bal_flag"] == 2) * (data["l_o3"] > 0) 
		print np.sum(select2)

		if np.sum(select_null*select2) > 0:

			if usebins[i] == None:
				n,bins_bals,patches=hist(data_to_plot[select_null*select2], facecolor=colors[1], alpha=0.5, log=LOG, label="BALs", normed=NORM, stacked=True)
			else:
				n,bins_bals,patches=hist(data_to_plot[select_null*select2], bins = usebins[i], facecolor=colors[1], alpha=0.5, log=LOG, label="BALs", normed=NORM, stacked=True)
			
			hist(data_to_plot[select_null*select1],bins=bins_bals, facecolor=colors[0], alpha=0.5, log=LOG, label="non-BALs", normed=NORM, stacked=True)
		title(data_labels[to_plot[i]], fontsize=18)
		if i == 0: float_legend(loc=2)

	savename = "hist_all.png"
	if select_z: savename = "zhist_all.png"

	savefig(savename, dpi=300)
	clf()


def make_scatter():
	figure(figsize=(12,18))
	scatter_to_plot = [("mbh_hbeta", "mbh_mg2"), ("mbh_hbeta", "l_o3"),
	("L5100", "ew_o3"), ("L5100", "l_o3"), ("L5100", "mbh_hbeta"),
	("ew_o3","l_o3")]

	for i in range(len(scatter_to_plot)):

		subplot(3,2,i+1)

		x = data[scatter_to_plot[i][0]]
		y = data[scatter_to_plot[i][1]]
		xerr = data["e_"+scatter_to_plot[i][0]]
		yerr = data["e_"+scatter_to_plot[i][1]]

		select_null = (x > 0) * (y > 0) 
		if select_mass: select_null *= (data["mbh_hbeta"] > 8.5) * (data["mbh_hbeta"] < 9.5)
		select1 = (data["bal_flag"] == 0) 
		select2 = (data["bal_flag"] >= 1) 

		scatter(x[select_null*select1], y[select_null*select1], c=colors[0])
		scatter(x[select_null*select2], y[select_null*select2], c=colors[1], s=100)
		errorbar(x[select_null*select2], y[select_null*select2], xerr=xerr[select_null*select2], yerr=yerr[select_null*select2], fmt='.', c=colors[1], linewidth=2)

		xlabel(data_labels[scatter_to_plot[i][0]], fontsize=18)
		ylabel(data_labels[scatter_to_plot[i][1]], fontsize=18)


	savefig("scatter_all.png", dpi=300)
	clf()

def make_hist_dem():

	figure(figsize=(12,18))

	to_plot2 = np.array(["radio_flag", "special_flag", "ew_o3", "ew_mg2", "L3000", "L5100"])
	logs = [False, False, True, True, False, False]
	usebins = [[0.5,1,1.5,2,2.5],np.arange(-0.5,4.5,0.5),np.arange(-1,3,0.1), np.arange(-1,3,0.1), np.arange(44,47,0.2), np.arange(44,47,0.2)]
	select1 = (data["bal_flag"] == 0) 
	select2 = (data["bal_flag"] >= 1) 

	for i in range(2):

		subplot(3,2,i+1)
		if i == 0: 
			select_null = (data[to_plot2[i]] > 0) 
		else:
			select_null = (data[to_plot2[i]] > -1) 

		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		hist(data[to_plot2[i]][select_null*select2], facecolor=colors[1], alpha=0.5, log=LOG, label="BALs", normed=NORM, bins=usebins[i])
		hist(data[to_plot2[i]][select_null*select1], facecolor=colors[0], alpha=0.5, log=LOG, label="non-BALs", normed=NORM, bins=usebins[i])	
		
		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

	select_dem = (data["special_flag"] > 0) 
	select_nodem = (data["special_flag"] == 0) 

	for i in range(2,4):

		print to_plot2[i]
		subplot(3,2,i+1)

		select_null = (data[to_plot2[i]] > 0) 
		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		if logs[i]: data_to_plot = np.log10(data[to_plot2[i]])

		print np.sum(select_null*select2*select_dem), np.sum(select_null*select2*select_nodem)


		hist(data_to_plot [select_null*select2*select_dem], facecolor=colors[1], alpha=0.5, log=LOG, label="BALs, DEM", normed=NORM, bins=usebins[i])
		hist(data_to_plot [select_null*select2*select_nodem], facecolor=colors[2], alpha=0.5, log=LOG, label="BALs, non-DEM", normed=NORM, bins=usebins[i])
		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

		subplot(3,2,i+3)
		hist(data_to_plot [select_null*select1*select_dem], facecolor=colors[0], alpha=0.5, log=LOG, label="non-BALs, DEM", normed=NORM, bins=usebins[i])
		hist(data_to_plot [select_null*select1*select_nodem], facecolor=colors[4], alpha=0.5, log=LOG, label="non-BALs, non-DEM", normed=NORM, bins=usebins[i])
		
		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

	savefig("hist_dem.png")


def make_hist_dem_nonbal():

	figure(figsize=(12,18))

	to_plot2 = np.array(["ew_o3", "ew_mg2", "L3000", "L5100", "l_o3", "l_o3", "mbh", "Lbol"])
	logs = [True, True, False, False, False, False, False, False]
	usebins = [np.arange(-1,3,0.1), np.arange(-1,3,0.1), np.arange(44,47,0.2), np.arange(44,47,0.2), np.arange(40,44,0.2), np.arange(40,44,0.2), np.arange(7.5,10,0.2), np.arange(44,47,0.2)]

	select_dem = (data["special_flag"] > 0) 
	select_nodem = (data["special_flag"] == 0) 

	for i in range(len(to_plot)):

		print to_plot2[i]
		subplot(4,2,i+1)

		select_null = (data[to_plot2[i]] > 0) 
		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		uselog = LOG

		data_to_plot = data[to_plot2[i]]
		if logs[i]: 
			data_to_plot = np.log10(data[to_plot2[i]])



		if np.sum(select_null*select_dem) == 0:
			uselog=False
		if np.sum(select_null*select_nodem) == 0:
			uselog=False
		try:
			hist(data_to_plot [select_null*select_dem], facecolor=colors[0], alpha=0.5, log=uselog, label="DEM", normed=NORM, bins=usebins[i])
			hist(data_to_plot [select_null*select_nodem], facecolor=colors[4], alpha=0.5, log=uselog, label="non-DEM", normed=NORM, bins=usebins[i])
		except:
			print "BA"


		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

	savefig("hist_dem_nonbal.png")


def make_hist_radio_nonbal():

	figure(figsize=(12,18))

	to_plot2 = np.array(["ew_o3", "ew_mg2", "L3000", "L5100", "l_o3", "l_o3", "mbh", "Lbol"])
	logs = [True, True, False, False, False, False, False, False]
	usebins = [np.arange(-1,3,0.1), np.arange(-1,3,0.1), np.arange(44,47,0.2), np.arange(44,47,0.2), np.arange(40,44,0.2), np.arange(40,44,0.2), np.arange(7.5,10,0.2), np.arange(44,47,0.2)]

	select_core = (data["radio_flag"] == 1)  
	select_lobe = (data["radio_flag"] == 2) 

	for i in range(len(to_plot)):

		print to_plot2[i]
		subplot(4,2,i+1)

		select_null = (data[to_plot2[i]] > 0) 
		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		uselog = LOG

		data_to_plot = data[to_plot2[i]]
		if logs[i]: 
			data_to_plot = np.log10(data[to_plot2[i]])



		if np.sum(select_null*select_core) == 0:
			uselog=False
		if np.sum(select_null*select_lobe) == 0:
			uselog=False
		try:
			hist(data_to_plot [select_null*select_core], facecolor=colors[0], alpha=0.5, log=uselog, label="Core", normed=NORM, bins=usebins[i])
			hist(data_to_plot [select_null*select_lobe], facecolor=colors[4], alpha=0.5, log=uselog, label="Lobe", normed=NORM, bins=usebins[i])
		except:
			print "BA"


		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

	savefig("hist_dem_radio.png")


def make_hist_radio():

	figure(figsize=(12,18))

	to_plot2 = np.array(["radio_flag", "special_flag", "ew_o3", "ew_mg2", "L3000", "L5100"])
	logs = [False, False, True, True, False, False]
	usebins = [[0.5,1,1.5,2,2.5],np.arange(-0.5,4.5,0.5),np.arange(-1,3,0.1), np.arange(-1,3,0.1), np.arange(44,47,0.2), np.arange(44,47,0.2)]
	select1 = (data["bal_flag"] == 0) 
	select2 = (data["bal_flag"] >= 1) 

	for i in range(2):

		subplot(3,2,i+1)
		if i == 0: 
			select_null = (data[to_plot2[i]] > 0) 
		else:
			select_null = (data[to_plot2[i]] > -1) 

		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		hist(data[to_plot2[i]][select_null*select2], facecolor=colors[1], alpha=0.5, log=LOG, label="BALs", normed=NORM, bins=usebins[i])
		hist(data[to_plot2[i]][select_null*select1], facecolor=colors[0], alpha=0.5, log=LOG, label="non-BALs", normed=NORM, bins=usebins[i])	
		
		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

	select_core = (data["radio_flag"] == 1)  
	select_lobe = (data["radio_flag"] == 2) 

	for i in range(2,4):

		print to_plot2[i]
		subplot(3,2,i+1)

		select_null = (data[to_plot2[i]] > 0) 

		if select_z: select_null *= (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1])

		if logs[i]: data_to_plot = np.log10(data[to_plot2[i]])

		hist(data_to_plot [select_null*select2*select_core], facecolor=colors[1], alpha=0.5, log=LOG, label="BALs, Core", normed=NORM, bins=usebins[i])
		hist(data_to_plot [select_null*select2*select_lobe], facecolor=colors[2], alpha=0.5, log=LOG, label="BALs, Lobe", normed=NORM, bins=usebins[i])
		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

		subplot(3,2,i+3)
		hist(data_to_plot [select_null*select1*select_core], facecolor=colors[0], alpha=0.5, log=LOG, label="non-BALs, Core", normed=NORM, bins=usebins[i])
		hist(data_to_plot [select_null*select1*select_lobe], facecolor=colors[4], alpha=0.5, log=LOG, label="non-BALs, Lobe", normed=NORM, bins=usebins[i])
		
		title(data_labels[to_plot2[i]], fontsize=18)
		float_legend()

	savefig("hist_radio.png")


make_hist_radio_nonbal()
# make_hist_dem_nonbal()
# #make_scatter()
# #make_hist()
# make_hist_dem()

#make_hist_radio()