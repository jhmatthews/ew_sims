#!/usr/bin/env python
import pylab as pl
import os, sys 
from constants import *
import numpy as np
import pyfits
import sdss_sub_data as sub 
from scipy.integrate import trapz
import py_plot_util as util



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

		elif datatype == "bat":

			self.seyfert = (data["TYPE"] == "Sy1") + (data["TYPE"] == "Sy1.2") + (data["TYPE"] == "Sy1.5") 
			self.sy1 = (data["TYPE"] == "Sy1")
			self.qso = (data["TYPE"] == "Quasar")
			self.z = (data["REDSHIFT"] >= 0)


def get_sdss_hst_matches(data, d_hst):
	
	plot_flags = np.zeros(len(data["z"]), dtype="int")
	plot_flags_for_hst = np.zeros(len(d_hst ["ra"]), dtype="int")
	map_for_hst = np.zeros(len(d_hst ["ra"]), dtype="int")
	map_for_sdss = np.zeros(len(data["z"]), dtype="int")

	for i in range(len(d_hst["ra"])):

		dra = np.fabs(d_hst["ra"][i] - data["ra"])
		ddec = np.fabs(d_hst["dec"][i] - data["dec"])

		dist = np.sqrt(dra**2 + ddec**2)

		#print dist[dist.argmin()]
		dz = np.fabs(data["z"][dist.argmin()] - d_hst["z"][i])

		dist_arcminute = 0.001
		#dist_arcminute = 1.0/3600.0
		if dist[dist.argmin()] < dist_arcminute and dz < 0.01:
			#print i, dist[dist.argmin()], dist.argmin(), data["z"][dist.argmin()], d_hst["SNR"][i]
			plot_flags[dist.argmin()]  = 1
			plot_flags_for_hst[i] = 1
			map_for_hst[i] = dist.argmin()
			map_for_sdss[dist.argmin()] = i


	plot_flags = plot_flags.astype(bool)

	return plot_flags, plot_flags_for_hst.astype(bool), map_for_sdss, map_for_hst


def BALnicity ( w, f, wline, norm_method = "SP", point=None):

	'''
	calculate the BALnicity index of a line
	
	:INPUT
		nu_line		float
					The frequency of the line at line centre
		
		nu_array		float
					array of frequencies
		
		spec_array	float
					array of intensities/fluxes
		
	:OUTPUT
		BALnicity index		float
		
	'''
	
	# we integrate over positive velocities as looking for blueshifted BALs

	# first we need to normalise to the continuum
	if norm_method == "fit":
		return -999
	elif norm_method == "SP":
		if point == None: 
			point = wline + 100.0
		norm = util.get_flux_at_wavelength( w, f, point)

	elif norm_method == "mean":
		p1 = util.get_flux_at_wavelength( w, f, point[0])
		p2 = util.get_flux_at_wavelength( w, f, point[0])
		norm = np.mean([p1, p2])
	

	fnorm = f / norm

	fv = 1 - (fnorm / 0.9)

	vel = C * (wline - w) / wline * 1e-5 	# velocity in KM/S

	constant = np.zeros(len(vel))
	for i, v in enumerate(vel):

		if v > 2000.0:
			if fnorm[i] < 0.9:
				select = (vel >= (v - 2000.0)) * (vel <= v)
				if fnorm[select].all() <= 0.9:
					constant[i] = 1.0

	integration_limits = (vel < 25000.) * (vel > 3000.)

	if np.sum(constant) == 0:
		return 0.0, norm

	BI = trapz(constant*fv[integration_limits], vel[integration_limits])	
	return BI, norm



# get the data
data = sub.get_sample("../data/catalog.dat")
d_hst = sub.get_hst("../data/hst_catalog.dat")

matches, hst_matches, map_for_sdss, map_for_hst = get_sdss_hst_matches(data, d_hst)

folders = np.loadtxt("../data/hst_catalog.dat", usecols=(1,), unpack=True, dtype="string")
hst_redshifts = d_hst["z"][hst_matches]
folder_matches = folders[hst_matches]
hst_snr = d_hst["SNR"][hst_matches]

bal_catalog = open("bal_hst_catalog.dat", "w")

lines = [1240.0, 1400.0, 1550.0]
points = [1350.0, 1350.0, 1700.0]
points = [ [1450.0, 1650.0] ]
lines = [1550.0]

bal_flags = np.zeros(len(folder_matches))
BIs = np.zeros( [len(folder_matches), len(lines)] )

bal_catalog.write("# Name Redshift RA DEC\n")

match_with_coadd = 0

c4_flags = np.zeros(len(folder_matches))

bal_names = np.loadtxt("BALs", dtype="string")

for i, f in enumerate(folder_matches):

#for i in range(100,200):

	f = folder_matches[i]

	print i

	# get the redshift from the quasar HST catalog
	z = hst_redshifts[i] 

	# find all the fitsfiles for this object
	fitsfiles_only = [fff for fff in os.listdir("hst_datapile/%s/" % f) if ".fits" in fff and "coadd" in fff and "all" in fff]
	fitsfiles = ["hst_datapile/%s/%s" % (f, fff) for fff in fitsfiles_only]

	plot = True

	if len(fitsfiles) > 0:
		match_with_coadd += 1

	# go through each fits file
	for j, fitsfile in enumerate(fitsfiles):

		# get the data
		d = pyfits.getdata(fitsfile)

		try:
			wrest = d["WAVE"] 
			w = wrest / (1. + z)
			flux = d["FLUX"]
			spectrum_wave=True

		except KeyError:
			spectrum_wave = False

		# try:
		# 	warray = d["WAVELENGTH"]
		# 	#print "!! ", fitsfiles[j], warray.shape
		# 	spectrum_wl = True
		# 	if len(warray) > 1:
		# 		wrest = d["WAVELENGTH"][1] 
		# 		w = wrest / (1. + z)
		# 		flux = d["FLUX"][1]
		# 	elif len(warray) == 1:
		# 		wrest = d["WAVELENGTH"][0]
		# 		w = wrest / (1. + z)
		# 		flux = d["FLUX"][0]
		# 	else:
		# 		spectrum_wl = False
		# except KeyError or IndexError:
		# 	spectrum_wl = False


		spectrum = spectrum_wave

		if spectrum: 

			ly_alpha_mask = (np.fabs(wrest - 1215.0) > 2.0) 
			flux_mask = (flux > 0.0) 
			#line_mask = (np.fabs(w - 1550.) < 300.0)

			mask = flux_mask * ly_alpha_mask
			wuse = w[mask]
			fuse = util.smooth(flux[mask])

			try:
				wmin = min(wuse)
				wmax = max(wuse)
			except ValueError:
				spectrum = False


		#spectrum = False


		if spectrum: 

			# flags to decide if we print out or not
			bal_flags_this_file = np.zeros(len(lines))

			BI = 0 

			if np.max(wuse) > 1600.0 and np.min(wuse) < 1500.0:
				c4_flags[i] = 1


			# for iline, line in enumerate(lines):


			# 	# are our lines and normalisation points in the spctrum?
			# 	do = (points[iline][0] > wmin) * (points[iline][1] < wmax)

			# 	if do:

			# 		#try:
			# 		BI, norm = BALnicity ( wuse, fuse, line, norm_method ="mean", point=points[iline])
			# 		#except:
			# 		#	print "failed"
			# 		norm_use = norm

			# 		if BI > 0:
			# 			bal_flags[i, iline] = 1
			# 			bal_flags_this_file[iline] = 1
			# 			BIs[i, iline] = BI
			# 			norm_use = norm
			# 			print BI


			#if plot==True:
			if False:
				pl.figure()

				lims = True

				try:
					ymax = np.max(fuse[(wuse > points[0][0]) * (wuse < points[0][1])])
				except:
					lims = False


				print type(wuse), type(fuse)
				pl.plot(wuse, fuse)	
				#pl.hlines([0.9*norm_use],1200,1700, linewidth=2)
				#pl.vlines( [1550.0 - (3e8 / C)*1550.0,1550.0 - (2.5e9 / C)*1550.0], 0,ymax, linewidth=2)

				#if lims:
				#	pl.ylim(0,ymax)
				pl.title("C4 flag is %i SNR is %8.4e" % (c4_flags[i], d_hst["SNR"][hst_matches][i]))

				#title("BIs: %.2f %.2f %.2f" % (BIs[i, 0], BIs[i, 1], BIs[i, 2]) )
				#pl.title("BIs = %.2f" % (BIs[i, 0]) )
				pl.xlim(600,1700)
				pl.semilogy()
				if c4_flags[i] == 0:
					folder_to_save = "cos_spectra_plots/no_c4/"
				else:
					folder_to_save = "cos_spectra_plots/with_c4/"

				pl.xlabel("Rest Wavelength (A)")
				pl.ylabel("Flux")
				pl.savefig("%s/%s_%i.png" % (folder_to_save, f, j))
				pl.clf()

				#sys.exit()

				print f, fitsfiles_only[j]

				plot = False

	select_bal_names = (f == bal_names)

	if np.sum(select_bal_names == 1):
		print f

		bal_flags[i] = 1

bal_catalog.close()

def get_dist(z):

	'''approximate distance in pc from redshift'''

	return (z * 3e5 / 70.0) * 1e6


from pylab import *


# things to plot from SDSS
plot_flags = matches

with_c4_plot_flags = map_for_hst[c4_flags.astype(bool)]
bal_plot_flags = map_for_hst[bal_flags.astype(bool)]


print match_with_coadd

d = get_dist(data["z"]) * PARSEC

select = selection(data)

times = [20.1, 201.0, 2010.]
fluxes = [-11., -12., -13.]

# calculate monochromatic fluxes
f1350 = np.log10(10.0**data["L1350"] / 4.0 / PI / d / d)
f3000 = np.log10(10.0**data["L3000"] / 4.0 / PI / d / d)
f5100 = np.log10(10.0**data["L5100"] / 4.0 / PI / d / d)

figure()
y1 = -14
y2 = -10

scatter(data["z"][select.z], f5100[select.z], label="SDSS DR7 Quasars",marker=".", edgecolors="None", facecolors="k", alpha=0.5)
scatter(data["z"][select.z*plot_flags], f5100[select.z*plot_flags], label="Also in COS archive",marker="o", facecolors="r", edgecolors="None", s=20)
scatter(data["z"][with_c4_plot_flags], f5100[with_c4_plot_flags], label="COS with CIV",marker="o", facecolors="r", edgecolors="None", s=20)


#semilogy()
ylim(y1, y2)
xlim(0,0.5)
xlabel("Redshift", fontsize=20)
ylabel(r"$F_{\lambda}$", fontsize=20)
#pretty_legend()
#vlines()

hlines(fluxes, 0, 1, linewidth=2)


#vlines(redshift_lims, -14,-10, linewidth=2)
#vlines(ly_alpha_lims, -14,-10, linewidth=2)
for i, t in enumerate(times):
	text(0.1, fluxes[i]+0.04, "%is" % t, fontsize=16)



savefig("fluxes_cos_feasibility_c4.png", dpi=300)
clf()

figure(figsize=(10,10))

subplot(2,2,1)
scatter(data["z"][select.z], data["edd_frac"][select.z], label="SDSS DR7 Quasars",marker=".", edgecolors="None", facecolors="k", alpha=0.5)
scatter(data["z"][select.z*plot_flags], data["edd_frac"][select.z*plot_flags], label="Also in COS archive",marker="o", facecolors="r", edgecolors="None", s=20)
scatter(data["z"][with_c4_plot_flags], data["edd_frac"][with_c4_plot_flags], label="Also in COS archive",marker="o", facecolors="b", edgecolors="None", s=20)
xlabel("Redshift", fontsize=20)
ylabel(r"$L/L_{edd}$", fontsize=20)
ylim(-3,2)
xlim(0,1)

subplot(2,2,2)
scatter(data["z"][select.z], data["mbh"][select.z], label="SDSS DR7 Quasars",marker=".", edgecolors="None", facecolors="k", alpha=0.5)
scatter(data["z"][select.z*plot_flags], data["mbh"][select.z*plot_flags], label="Also in COS archive",marker="o", facecolors="r", edgecolors="None", s=20)
scatter(data["z"][with_c4_plot_flags], data["mbh"][with_c4_plot_flags], label="COS with CIV",marker="o", facecolors="b", edgecolors="None", s=20)
xlabel("Redshift", fontsize=20)
ylabel(r"$M_{BH}$", fontsize=20)
ylim(7,11)
xlim(0,1)



subplot(2,2,3)
scatter(data["edd_frac"][select.z], data["mbh"][select.z], label="SDSS DR7 Quasars",marker=".", edgecolors="None", facecolors="k", alpha=0.5)
scatter(data["edd_frac"][select.z*plot_flags], data["mbh"][select.z*plot_flags], label="Also in COS archive",marker="o", facecolors="r", edgecolors="None", s=20)
scatter(data["edd_frac"][with_c4_plot_flags], data["mbh"][with_c4_plot_flags], label="COS with CIV",marker="o", facecolors="b", edgecolors="None", s=20)
xlabel(r"$L/L_{edd}$", fontsize=20)
ylabel(r"$M_{BH}$", fontsize=20)
ylim(7,11)
xlim(-3,2)

subplot(2,2,4)
scatter(d_hst["z"], d_hst["SNR"], label="Also in COS archive",marker="o", facecolors="r", edgecolors="None", s=20)
scatter(d_hst["z"][hst_matches][c4_flags.astype(bool)], d_hst["SNR"][hst_matches][c4_flags.astype(bool)], label="COS with CIV",marker="o", facecolors="b", edgecolors="None", s=20)

xlabel("Redshift", fontsize=20)
ylabel(r"Median $SNR$", fontsize=20)
ylim(0,100)
xlim(0,1)
savefig("v_redshift_c4.png", dpi=300)

clf()


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
pretty.pretty_legend()
subplots_adjust(hspace=0, wspace=0, bottom=0.1,top=0.95)
savefig("bals_scatter.png", dpi=300)













