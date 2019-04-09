from scipy import *
from numpy import *
from numpy.fft import *
from math import *
from pylab import *
import scipy.signal
import scipy.stats
import os
import pdb
import random
import sys
import time
import pickle
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#######################################################
# bins data and calculates average/STD within bins
def binned_averages_std(stdpData,bins):
        binnedData = zeros((len(bins)-1,4))
        for i in range(len(bins)-1):
                mask = (stdpData[:,0] >= bins[i]) & (stdpData[:,0]< bins[i+1])
                binnedData[i,0] = average(stdpData[mask][:,0])
                binnedData[i,1] = average(stdpData[mask][:,1])
                binnedData[i,2] = std(stdpData[mask][:,1])
                binnedData[i,3] = std(stdpData[mask][:,0])
        return binnedData
        


#######################################################
# load experimental data
dataLocation = '/home/mgraupe/theobio/camkII/simplest_model/laurent_venance/experimental_data/'
stdp1Hz = np.loadtxt('STDP_1Hz_100pairings.dat')
stdp3Hz = np.loadtxt('STDP_2.5-3Hz_100pairings.dat')
stdp5Hz = np.loadtxt('STDP_5Hz_100pairings.dat')
stdp10Hz = np.loadtxt('STDP_10Hz_100pairings.dat')

# binning of the raw data
bins1 = array([-220,-190,-40,-30,-20,-10,0,15,20,30,40,250])
stdp1binned = binned_averages_std(stdp1Hz,bins1)

bins3 = array([-110,-80,-60,-40,-20,0,20,60,80,110])
stdp3binned = binned_averages_std(stdp3Hz,bins3)

bins5 = array([-110,-80,-50,-20,0,20,61])
stdp5binned = binned_averages_std(stdp5Hz,bins5)

bins10 = array([-81,-60,-40,0,20,40,71])
stdp10binned = binned_averages_std(stdp10Hz,bins10)

# saved binned data
np.savetxt('STDP_1Hz_100pairings_binned.dat',stdp1binned)
np.savetxt('STDP_2.5-3Hz_100pairings_binned.dat',stdp3binned)
np.savetxt('STDP_5Hz_100pairings_binned.dat',stdp5binned)
np.savetxt('STDP_10Hz_100pairings_binned.dat',stdp10binned)


#######################################################
# plot data
fig_width = 12 # width in inches
fig_height = 8  # height in inches
fig_size =  [fig_width,fig_height]
params = {'axes.labelsize': 14,
          'axes.titlesize': 13,
          'font.size': 11,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11,
          'figure.figsize': fig_size,
          'savefig.dpi' : 600,
          'axes.linewidth' : 1.3,
          'ytick.major.size' : 4,      # major tick size in points
          'xtick.major.size' : 4      # major tick size in points
          #'edgecolor' : None
          #'xtick.major.size' : 2,
          #'ytick.major.size' : 2,
          }
rcParams.update(params)

# set sans-serif font to Arial
rcParams['font.sans-serif'] = 'Arial'

# create figure instance
fig = plt.figure()


# define sub-panel grid and possibly width and height ratios
gs = gridspec.GridSpec(2, 2
                       #width_ratios=[1,1.2],
                       #height_ratios=[1,1]
                       )

# define vertical and horizontal spacing between panels
gs.update(wspace=0.2,hspace=0.3)

fig.suptitle('STDP data, Laurent Venance lab (100 regular-spaced spike-pairs)',fontsize=14)
# possibly change outer margins of the figure
#plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.18)

# sub-panel enumerations
#plt.figtext(0.06, 0.92, 'A',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.47, 0.92, 'B',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.06, 0.47, 'C',clip_on=False,color='black', weight='bold',size=22)
#plt.figtext(0.47, 0.47, 'D',clip_on=False,color='black', weight='bold',size=22)


# panel 0 #######################################################
ax0 = plt.subplot(gs[0])

# title
ax0.set_title('1 Hz')

# diplay of data
ax0.axhline(y=100,ls='--',color='0.7',lw=2)
ax0.axvline(x=0,ls='--',color='0.7',lw=2)
ax0.plot(stdp1Hz[:,0],stdp1Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
ax0.errorbar(stdp1binned[:,0],stdp1binned[:,1],yerr=stdp1binned[:,2],fmt='o-',markeredgecolor='green')

# removes upper and right axes 
# and moves left and bottom axes away
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')

ax0.set_xlim(-250,250)
# legends and labels
#plt.legend(loc=1,frameon=False)

#plt.xlabel(r'$\Delta t$ (ms)')
plt.ylabel('change in synaptic strength')


ax_inset = plt.axes((0.35, 0.71, 0.15, 0.2))
#plt.hist(residuals, fc='0.8',ec='w', lw=2)
ax_inset.axhline(y=100,ls='--',color='0.7',lw=2)
ax_inset.axvline(x=0,ls='--',color='0.7',lw=2)
ax_inset.plot(stdp1Hz[:,0],stdp1Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
ax_inset.errorbar(stdp1binned[:,0],stdp1binned[:,1],yerr=stdp1binned[:,2],fmt='o-',markeredgecolor='green')
ax_inset.set_xlim(-40,40)

ax_inset.spines['top'].set_visible(False)
ax_inset.spines['right'].set_visible(False)
ax_inset.spines['bottom'].set_position(('outward', 10))
ax_inset.spines['left'].set_position(('outward', 10))
ax_inset.yaxis.set_ticks_position('left')
ax_inset.xaxis.set_ticks_position('bottom')

ax_inset.yaxis.set_major_locator(MaxNLocator(3))
ax_inset.xaxis.set_major_locator(MaxNLocator(3))
#plt.xticks([-0.5, 0, 0.5],[-0.5, 0, 0.5], size=6)
#plt.xlim((-0.5, 0.5))
#plt.yticks([5, 10], size=6)
#plt.xlabel("residuals", size=8)
#ax_inset.xaxis.set_ticks_position("none")
#ax_inset.yaxis.set_ticks_position("left")


# panel 1 #############################################
ax1 = plt.subplot(gs[1])

# title
ax1.set_title('2.5 - 3 Hz')

# diplay of data
ax1.axhline(y=100,ls='--',color='0.7',lw=2)
ax1.axvline(x=0,ls='--',color='0.7',lw=2)
ax1.plot(stdp3Hz[:,0],stdp3Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
ax1.errorbar(stdp3binned[:,0],stdp3binned[:,1],yerr=stdp3binned[:,2],fmt='o-',markeredgecolor='green')

# removes upper and right axes 
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

ax1.set_xlim(-110,110)
# legends and labels
#plt.legend(loc=1,frameon=False)

#plt.xlabel(r'$\Delta t$ (ms)')
#plt.ylabel('change in synaptic strength')


# third sub-plot #######################################################
ax2 = plt.subplot(gs[2])


# title
ax2.set_title('5 Hz')

# diplay of data
ax2.axhline(y=100,ls='--',color='0.7',lw=2)
ax2.axvline(x=0,ls='--',color='0.7',lw=2)
ax2.plot(stdp5Hz[:,0],stdp5Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
ax2.errorbar(stdp5binned[:,0],stdp5binned[:,1],yerr=stdp5binned[:,2],fmt='o-',markeredgecolor='green')

# removes upper and right axes 
# and moves left and bottom axes away
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_position(('outward', 10))
ax2.spines['left'].set_position(('outward', 10))
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')

ax2.set_xlim(-110,110)
# legends and labels
#plt.legend(loc=1,frameon=False)

plt.xlabel(r'$\Delta t$ (ms)')
plt.ylabel('change in synaptic strength')

# fourth sub-plot #######################################################
ax3 = plt.subplot(gs[3])


# title
ax3.set_title('10 Hz')

# diplay of data
ax3.axhline(y=100,ls='--',color='0.7',lw=2)
ax3.axvline(x=0,ls='--',color='0.7',lw=2)
ax3.plot(stdp10Hz[:,0],stdp10Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5',clip_on=False)
ax3.errorbar(stdp10binned[:,0],stdp10binned[:,1],yerr=stdp10binned[:,2],fmt='o-',markeredgecolor='green')

# removes upper and right axes 
# and moves left and bottom axes away
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_position(('outward', 10))
ax3.spines['left'].set_position(('outward', 10))
ax3.yaxis.set_ticks_position('left')
ax3.xaxis.set_ticks_position('bottom')

ax3.set_xlim(-110,110)
# legends and labels
#plt.legend(loc=1,frameon=False)

plt.xlabel(r'$\Delta t$ (ms)')
#plt.ylabel('change in synaptic strength')

## save figure ############################################################
fname = os.path.basename(__file__)[:-3]

savefig(fname+'.png')
savefig(fname+'.pdf')