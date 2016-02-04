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



def plot_contour(sim):

	'''
	sim 
	instance of class simulation 
	'''

	cm_use=get_viridis()
	# below here it's just plotting...
	figure(figsize=(13,7))
	subplot(1,2,1)

	sim.ks_p_value = np.ma.masked_array(sim.ks_p_value, mask=(sim.ks_p_value == 0) )
	sim.f_bal = np.ma.masked_array(sim.f_bal, mask=(sim.f_bal == 0) )

	contourf(sim.thetamax, sim.thetamin, np.log10(sim.ks_p_value), extend='both', levels=np.arange(-5,0,0.1), cmap=cm_use)
	#contour(thetamax, thetamins, np.log10(ks_test), levels=np.log10(np.array([0.001,0.05,0.32])), c="k")
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title('''$\log~(p$-value$)$. If the p-value is high, 
	then we cannot reject the hypothesis that the distributions 
	of the two samples are the same''', fontsize=14)

	subplot(1,2,2)
	contourf(sim.thetamax, sim.thetamin, sim.f_bal, extend='both', levels=np.arange(0,1,0.02), cmap=cm_use)
	#contour(thetamax, thetamins, f_bal)
	colorbar()
	ylabel(r"$\theta_{min}$", fontsize=20)
	xlabel(r"$\theta_{max}$", fontsize=20)
	title("$f_{BAL}$", fontsize=20)

	subplots_adjust(left=0.1,right=0.97,top=0.80)
	savefig("contour.png", dpi=200)

	return 0


def plot_contour2(sim, select):

	'''
	sim 
	instance of class simulation 
	'''

	d_use = sim.data["ew_o3"][select.general*select.mgbal]

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
	savefig("contour2.png", dpi=200)

	return 0
