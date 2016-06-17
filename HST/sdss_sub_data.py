#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np

# column names used in the dictionary used to access the SDSS data. Should match with the column map
column_names = np.array(["z","Lbol", "bal_flag", "radio_flag", 
"L5100", "L3000", "L1350", "l_o3", 
"ew_o3", "l_mg2", "ew_mg2", 
"l_c4", "ew_c4", "mbh_hbeta", "mbh_mg2", "mbh", "edd_frac","special_flag", "ra", "dec", "rflux","ew_fe_ha","ew_fe_hb", "ew_hb", "fwhm_hb"])


label_columns = np.array(["$z$","$L_{bol}$", "bal flag", "radio flag", 
"$L_{5100}$", "$L_{3000}$", "$L_{1350}$", r"L[O~\textsc{iii}]", 
r"EW[O~\textsc{iii}]", r"$L$[Mg~\textsc{ii}]", r"EW[Mg~\textsc{ii}]", 
r"$L$[C~\textsc{iv}]", r"EW[C~\textsc{iv}]", r"$M_{BH}(H\beta)$", r"$M_{BH}({\rm Mg~textsc{ii}})$",
"$M_{BH}$","$L/L_{edd}$","special flag", "RA", "Dec", "rflux",r"EW[Fe~\textsc{ii}]",r"EW[Fe~\textsc{ii}]"
, r"FWHM[H$\beta$]", r"EW[H$\beta$]"])


column_map = np.array([3, 11, 13, 14, 18, 20, 22, 72, 74, 85, 87, 105, 107, 126, 130, 138, 140, 141, 1, 2, 15, 48, 76, 59, 57])
column_map_with_errors = np.zeros(2*len(column_map), dtype=int)
for i in range(0,len(column_map_with_errors),2):

	icol = i / 2
	column_map_with_errors[i] = column_map[icol]
	column_map_with_errors[i+1] = column_map[icol]+1


def get_sample(fname, errors=True):

	'''
	Read the data and place in a dictionary object
	'''

	colmap = column_map_with_errors

	data = np.loadtxt(fname, unpack=True, usecols=colmap)
	data_dict = dict()


	for i in range(0,len(colmap),2):

		data_dict[column_names[i/2]] = data[i]
		data_dict["e_"+column_names[i/2]] = data[i+1]

	return data_dict


def make_label_map():

	'''
	map the labels to the column names for plotting and so on 
	'''

	labels = dict()

	for i in range(len(column_names)):
		labels[column_names[i]] = label_columns[i]

	return labels


def get_hst(fname):

	'''
	Read the data and place in a dictionary object
	'''

	colmap = ["ra", "dec", "z", "SNR"]
	#datatype = ["float", "float"]


	data_dict = dict()

	data = np.loadtxt(fname, unpack=True, usecols=(2,3,8, 10), dtype="string")

	data[2][(data[2] == "..")] = "-999"
	data[2][(data[2] == "..")] = "-999"




	for i in range(len(colmap)):

		#if datatype[i] == "float":
		data_dict[colmap[i]] = data[i].astype(float)
		#else:
		#	data_dict[colmap[i]] = data[i]

	return data_dict