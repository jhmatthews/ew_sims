#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import pretty 
import mynormalize
import sdss_sub_data as sub 
# from colorbar import DraggableColorbar as drbar

cc = pretty.get_colors()


class selection:

	'''
	a class which contains a number of Boolean arrays used for easy
	selection of subsamples in the SDSS data 
	'''

	def __init__(self, data, datatype="sdss"):

		if datatype == "sdss":
			self.nonbal = (data["bal_flag"] == 0) 
			self.mgbal = (data["bal_flag"] == 2)
			#self.z = (data["z"] > redshift_lims[0]) * (data["z"] < redshift_lims[1]) * (data["z"] < 2)
			self.z = (data["z"] < 1)
			self.has_o3 = (data["ew_o3"] > 0)
			self.mass = None
			self.general = self.has_o3 * self.z
			self.dem = (data["special_flag"] > 0)
			self.redshiftA = (data["z"] > 0.35) * (data["z"] < 0.83)
			self.redshiftB= (data["z"] > 1.45) * (data["z"] < 2.28)

		elif datatype == "bat":

			self.seyfert = (data["TYPE"] == "Sy1") + (data["TYPE"] == "Sy1.2") + (data["TYPE"] == "Sy1.5") 
			self.sy1 = (data["TYPE"] == "Sy1")
			self.qso = (data["TYPE"] == "Quasar")
			self.z = (data["REDSHIFT"] >= 0)


def figure_init():

	figure(figsize=(6.5,8.5))
	pretty.big_tick_labels(16)
	pretty.long_ticks()
	pretty.set_pretty()
	cc = pretty.get_colors()

	return 0

def make_subplot(xstring, ystring, data, select, loc=211, sample="A"):

	subplot(loc)

	if sample == "A":
		nonbals = select.redshiftA *select.nonbal
		bals = select.mgbal * select.redshiftA
		bal_label="LoBALs, S11"
		icolor = 2
		alpha = 1
	else:
		nonbals = select.redshiftB * select.nonbal
		bals = (data["bal_flag"] > 0) * select.redshiftB
		bal_label="BALs, S11"
		icolor = 0
		alpha = 0.5

	print "%i BALs in sample %s" % (np.sum(bals), sample)
	print "%s means are non-BAL: %8.4e BAL: %8.4e" % (xstring, np.mean(data[xstring][nonbals]), np.mean(data[xstring][bals]))
	print "%s means are non-BAL: %8.4e BAL: %8.4e" % (ystring, np.mean(data[ystring][nonbals]), np.mean(data[ystring][bals]))

	# edd_bins = np.arange(-3,2,0.05)
	# mbins = np.arange(7,11,0.05)
	#hist2d(data["edd_frac"], data["mbh"], bins=[edd_bins, mbins], normed=True, cmap="Greys")
	scatter(data[xstring][nonbals], data[ystring][nonbals], 
	        label="S11 Quasars",marker=".", edgecolors="None", 
	        facecolors="k", alpha=0.5)

	yerr = data["e_edd_frac"][bals]
	xerr = data["e_mbh"][bals]

	#errorbar(data["edd_frac"][mgbals], data["mbh"][mgbals],
	#xerr=xerr,yerr=yerr, fmt="o", ecolor=cc[2], markerfacecolor=cc[2], markeredgecolor="None", markersize=7, label="Mg BALs, S10")

	scatter(data[xstring][bals], data[ystring][bals],
	        c=cc[icolor], facecolor=cc[icolor], edgecolor="None",label=bal_label, alpha=alpha)

	xlim(-3,2)
	ylim(7,11)
	#ylabel(r"$\log [ M_{BH} ( M_\odot) ]$", fontsize=20)
	pretty.float_legend()
	text(-2.5,7.5, sample, fontsize=26)


# specify your variables here
XSTRING = "edd_frac"
YSTRING = "mbh"
xlims = (-3,2)
ylims = (7,11)


data = sub.get_sample("../data/catalog.dat")
select = selection(data)

figure_init()

make_subplot(XSTRING, YSTRING, data, select, loc=211, sample="A")

ylabel(r"$\log [ M_{BH} ( M_\odot) ]$", fontsize=20)

# bins = np.array([np.arange(-4,2.1,0.05), np.arange(7,11.1,0.05)])
# cmap = pretty.get_viridis()
# cmap="Greys"
#cnts, x, y, img = hist2d(data["edd_frac"][select.z], data["mbh"][select.z],bins=bins,cmap=cmap)
# yerr = data["e_edd_frac"][bal_plot_flags]
# xerr = data["e_mbh"][bal_plot_flags]
# # scatter(data["edd_frac"][bal_plot_flags], data["mbh"][bal_plot_flags], 
# # 	    label="HST Selected HiBALs",marker="o", edgecolors="None", 
# # 	    facecolors="r")
# errorbar(data["edd_frac"][bal_plot_flags], data["mbh"][bal_plot_flags],
# 	     xerr=xerr,yerr=yerr, fmt="o", ecolor=cc[1], markerfacecolor=cc[1], markeredgecolor="None", markersize=7, label="HiBALs, HST")


make_subplot(XSTRING, YSTRING, data, select, loc=212, sample="B")

# cbar = colorbar(orientation="horizontal")
# cbar.set_label("Number in S10 Catalog Bin", fontsize=16)
# cbar.set_norm(mynormalize.MyNormalize(vmin=0,vmax=1,stretch='linear'))
# #cbar =drbar(cbar,img)


xlabel(r"$\log [L_{bol} / L_{Edd}]$", fontsize=20)
ylabel(r"$\log [ M_{BH} ( M_\odot) ]$", fontsize=20)

# finalie and save figure
subplots_adjust(hspace=0.1, wspace=0.1, right=0.97, bottom=0.1,top=0.95)
savefig("%s_v_%s.png" % (XSTRING, YSTRING), dpi=300)