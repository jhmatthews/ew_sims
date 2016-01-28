#!/usr/bin/env python
from pylab import *
import py_read_output as r 
import py_plot_output as pl
import py_plot_util as util 
import os, sys 
from constants import *
import numpy as np
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = figure()
ax = fig.add_subplot(111, projection='3d')

n = 3e3

costhetas = np.power((np.random.random(size=n)),1./3.)
sinthetas = np.sqrt (1. - costhetas * costhetas)
thetas = np.arccos(costhetas)

phi = np.random.random(size=n) * 2.0 * np.pi


x =  np.cos(phi) * sinthetas
y =  np.sin(phi) * sinthetas
z =  costhetas



xlim(-1,1)
ylim(-1,1)
ax.set_zlim(-1,1)

ax.scatter(x, y, z)
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')

show()