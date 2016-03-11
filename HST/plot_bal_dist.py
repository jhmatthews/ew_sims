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
from colorbar import DraggableColorbar as drbar

figure()
pretty.big_tick_labels(16)
pretty.long_ticks()
pretty.set_pretty()
cc = pretty.get_colors()
bal_plot_flags = map_for_hst[bal_flags.astype(bool)]

scatter(data["edd_frac"][select.z], data["mbh"][select.z], 
	    label="S10 Quasars",marker=".", edgecolors="None", 
	    facecolors="k", alpha=0.5)


bins = np.array([np.arange(-4,2.1,0.05), np.arange(7,11.1,0.05)])
cmap = pretty.get_viridis()
cmap="Greys"
#cnts, x, y, img = hist2d(data["edd_frac"][select.z], data["mbh"][select.z],bins=bins,cmap=cmap)


yerr = data["e_edd_frac"][bal_plot_flags]
xerr = data["e_mbh"][bal_plot_flags]

# scatter(data["edd_frac"][bal_plot_flags], data["mbh"][bal_plot_flags], 
# 	    label="HST Selected HiBALs",marker="o", edgecolors="None", 
# 	    facecolors="r")

errorbar(data["edd_frac"][bal_plot_flags], data["mbh"][bal_plot_flags],
	     xerr=xerr,yerr=yerr, fmt="o", ecolor=cc[1], markerfacecolor=cc[1], markeredgecolor="None", markersize=7, label="HiBALs, HST")

xlim(-3,2)
ylim(7,11)
mgbals = select.general*select.mgbal
yerr = data["e_edd_frac"][mgbals]
xerr = data["e_mbh"][mgbals]

errorbar(data["edd_frac"][mgbals], data["mbh"][mgbals],
	     xerr=xerr,yerr=yerr, fmt="o", ecolor=cc[2], markerfacecolor=cc[2], markeredgecolor="None", markersize=7, label="Mg BALs, S10")

# cbar = colorbar(orientation="horizontal")
# cbar.set_label("Number in S10 Catalog Bin", fontsize=16)
# cbar.set_norm(mynormalize.MyNormalize(vmin=0,vmax=1,stretch='linear'))
# #cbar =drbar(cbar,img)


xlabel(r"$\log [L_{bol} / L_{Edd}]$", fontsize=20)
ylabel(r"$\log [ M_{BH} ( M_\odot) ]$", fontsize=20)
pretty.float_legend()
subplots_adjust(hspace=0, wspace=0, bottom=0.1,top=0.95)
