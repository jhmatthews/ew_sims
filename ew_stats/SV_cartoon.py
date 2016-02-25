'''
SV_cartoon.py 

This script creates the cartoon figure for 
the 2014 CV paper. The data were scraped from Dexter.
'''

import matplotlib.mlab as mlab
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pretty import *
from pylab import *

#set_pretty()
XKCD = False

if XKCD:
	plt.xkcd()
else:
	set_pretty()


#mpl.rcParams['pdf.use14corefonts'] = True
#mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#mpl.rc('text', usetex=True)

def draw_cartoon(thetamin, thetamax, ax, text=True):

	diskmax=30.0
	xmin=6.0
	xmax=25.0
	rstar=2.0

	xwindmax=43.0
	zwindmax=(xwindmax-xmin)*np.tan(np.radians(90.0-thetamin))


	posx=np.arange(xmin,xwindmax,(xmax-xmin)/10.0)
	ymin1=[]
	ymax1=[]

	x=np.array([xmin,xmax,xwindmax])
	y1=np.array([0.0,(xmax-xmin)*np.tan(np.radians(90.0-thetamin)),(xwindmax-xmin)*np.tan(np.radians(90.0-thetamin))])
	y2=np.array([0.0,0.0,(xwindmax-xmax)*np.tan(np.radians(90.0-thetamax))])

	lsize=14
	tsize=16 

	ax.plot(x,y1,x,y2,color='black')
	ax.set_xlim([-1.0*xwindmax,xwindmax])
	fill_between(x,y1, y2, facecolor='0.9', interpolate=True)
	plot([0.0,diskmax],[0.0,0.0],color='k',linewidth=3)
	plot(0.0,0.0,'o',color='w',ms=20.0)

	''' This is for the radiating sources'''
	#ax.annotate(r'$\rm{Central~BH}$', xy=(-0.5,0.6),xytext=(-17,8),fontsize=lsize,arrowprops=dict(arrowstyle="->"),)
	if text:
		ax.annotate(r'$\rm{Disc}$', xy=(5.0,0.1),xytext=(-5,3),fontsize=lsize,arrowprops=dict(arrowstyle="->"),)
		ax.annotate(r'$\rm{BH}$', xy=(-0.5,0.3),xytext=(-6,1.5),fontsize=lsize,arrowprops=dict(arrowstyle="->"),)
		ax.text(-2,8,"Symmetry Axis", rotation=90)
		ax.annotate('',xy=(29.5,2), xycoords='data',xytext=(20.6,3.5),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3 ,rad=-0.5"))
		ax.annotate('',xy=(16.2,2), xycoords='data',xytext=(0,4),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3 ,rad=-0.5"))


		'''This is for the wind parameters'''
		ax.plot([0,0],[0.8,10],'--',color='k')

	if text:
		if XKCD:
			ax.text(5,5,r"Type 1 Quasars", fontsize=12)
			ax.text(24,4,r"BAL Quasars'", fontsize=12)
			ax.text(28,0.7,r"`Obscured Quasars'", fontsize=12)
		else:
			ax.text(5,5,r"$0 \geq \theta<\theta_{min}$, `Type 1 Quasars'", fontsize=12)
			ax.text(24,4,r"$\theta_{min} \geq \theta<\theta_{max}$, `BAL Quasars'", fontsize=12)
			ax.text(28,0.7,r"$\theta>\theta_{max}$, `Obscured Quasars'", fontsize=12)
	else:
		#ax.text(-8,3,r"$\theta_{min} = %i^\circ, \theta_{max} = %i^\circ$" % (thetamin, thetamax), fontsize=16, rotation=90)
		ax.text(-8,7,r"$%i^\circ$ to $%i^\circ$" % (thetamin, thetamax), fontsize=16, rotation=90)

	ax.set_xlim(-3,45)
	ax.set_ylim(-1,11)
	axis('off')

fig=figure()
ax=fig.add_subplot(1,1,1)
savefig('fig2_cartoon.png',bbox_inches='tight', dpi=300)
#plt.savefig('../../figures/fig2_cartoon.eps',bbox_inches='tight')



