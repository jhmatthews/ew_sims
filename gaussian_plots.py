#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
from ew_stats_main import *
from SV_cartoon import draw_cartoon


thetamins = [30,45,60,75]
thetamaxes = [40,60,75,90]


figure(figsize=(8,10))
NPTS = 1000000
ew = np.random.normal( loc=10, scale=5, size=NPTS)
for i, thmin in enumerate(thetamins):

	thmax = thetamaxes[i]

	costhetas, bal_flags = get_mock_angles(thmin, NPTS, thmax)

	bins = np.arange(0,200,1)
	bins[0] = 0.0001
	bins = np.arange(-3,5,0.05)

	subplot(4,2,(2*i)+2)
	hist( np.log10(ew[(bal_flags == "q")] / costhetas[(bal_flags == "q")]), log=True, normed=True, bins=bins, label="non-BALs", alpha=0.6)
	hist( np.log10(ew[(bal_flags == "b")] / costhetas[(bal_flags == "b")]), log=True, normed=True, bins = bins, label="BALs", alpha=0.4)
	n, bin_edges = histogram( np.log10(ew), normed=True, bins=bins)

	plot(bins[:-1], n, c="r", label=r"$g({\rm EW})$", linewidth=1.5)

	ylabel("Frequency", fontsize=12)


	if i == 0: float_legend()

	xlim(0,3)


	subplot(4,2,(2*i)+1)
	draw_cartoon(thmin, thmax, gca(), text = False)


subplot(4,2,(2*i)+2)
xlabel(r"$\log [{\rm EW(\AA)}]$", fontsize=20)
subplots_adjust(left = 0.08, right = 0.95, top=0.95, bottom=0.1)

savefig("ew_with_cartoon.png", dpi=200)
clf()



