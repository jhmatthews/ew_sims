# #!/usr/bin/env python
# from pylab import *
# import py_read_output as r 
# import py_plot_output as pl
# import py_plot_util as util 
# import os, sys 
# from constants import *
# import numpy as np
# import pretty 
# import mynormalize
# from colorbar import DraggableColorbar as drbar

figure(figsize=(6.5,8.5))
pretty.big_tick_labels(16)
pretty.long_ticks()
pretty.set_pretty()
cc = pretty.get_colors()
#bal_plot_flags = map_for_hst[bal_flags.astype(bool)]

subplot(211)

edd_bins = np.arange(-3,2,0.05)
mbins = np.arange(7,11,0.05)

redshift = (data["z"] > 0.35) * (data["z"] < 0.83)

#hist2d(data["edd_frac"], data["mbh"], bins=[edd_bins, mbins], normed=True, cmap="Greys")
scatter(data["edd_frac"][redshift], data["mbh"][redshift], 
	    label="S11 Quasars",marker=".", edgecolors="None", 
	    facecolors="k", alpha=0.5)

mgbals = select.mgbal
yerr = data["e_edd_frac"][mgbals]
xerr = data["e_mbh"][mgbals]

#errorbar(data["edd_frac"][mgbals], data["mbh"][mgbals],
#	     xerr=xerr,yerr=yerr, fmt="o", ecolor=cc[2], markerfacecolor=cc[2], markeredgecolor="None", markersize=7, label="Mg BALs, S10")

scatter(data["edd_frac"][mgbals], data["mbh"][mgbals],
	     c=cc[2], facecolor=cc[2], edgecolor="None",label="LoBALs, S11", alpha=1)


xlim(-3,2)
ylim(7,11)
ylabel(r"$\log [ M_{BH} ( M_\odot) ]$", fontsize=20)
pretty.float_legend()
text(-2.5,7.5, "A", fontsize=26)

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


subplot(212)

redshift = (data["z"] > 1.45) * (data["z"] < 2.28)

#hist2d(data["edd_frac"], data["mbh"], bins=[edd_bins, mbins], normed=True, cmap="Greys")
scatter(data["edd_frac"][redshift*select.nonbal], data["mbh"][redshift*select.nonbal], 
	    label="S11 Quasars",marker=".", edgecolors="None", 
	    facecolors="k", alpha=0.5)

bals = (data["bal_flag"] == 1) * redshift
yerr = data["e_edd_frac"][bals]
xerr = data["e_mbh"][bals]
#scatter(data["edd_frac"][bals], data["mbh"][bals],
#	    c=cc[0])

scatter(data["edd_frac"][bals], data["mbh"][bals],
	     c=cc[0], facecolor=cc[0], edgecolor="None",label="BALs, S11", alpha=0.4)

xlim(-3,2)
ylim(7,11)
# cbar = colorbar(orientation="horizontal")
# cbar.set_label("Number in S10 Catalog Bin", fontsize=16)
# cbar.set_norm(mynormalize.MyNormalize(vmin=0,vmax=1,stretch='linear'))
# #cbar =drbar(cbar,img)

text(-2.5,7.5, "B", fontsize=26)
xlabel(r"$\log [L_{bol} / L_{Edd}]$", fontsize=20)
ylabel(r"$\log [ M_{BH} ( M_\odot) ]$", fontsize=20)
pretty.float_legend()
subplots_adjust(hspace=0.1, wspace=0.1, right=0.97, bottom=0.1,top=0.95)
savefig("bals_2x2_scatter.png", dpi=300)