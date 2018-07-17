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
from matplotlib import rcParams, font_manager
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import synapticChange
from timeAboveThreshold.timeAboveThreshold import *
from matplotlib.pyplot import gca

##########################################################
def calculateSnapticeChange(deltaT,frequency,interval,tat):

    interval = 1./frequency
    
    (timeD,timeP) = tat.spikePairFrequency(deltaT-synapticChange.D,frequency)
        
    # average potentiation and depression rates
    GammaP = synapticChange.gammaP*timeP/interval
    GammaD = synapticChange.gammaD*timeD/interval
    # rhoBar: average value of rho in the limit of a very long protocol equivalent to the minimum of the quadratic potentia
    rhoBar = GammaP/(GammaP + GammaD)
    # sigmaRhoSquared L standard deviation of rho in the same limit
    sigmaRhoSquared = (timeP/interval + timeD/interval)*(synapticChange.sigma**2)/(GammaP + GammaD)
    # tauEff : characteristic time scale of the temporal evolution of the pdf of rho
    tauEff = synapticChange.tau/(GammaP + GammaD)
    #
    # UP and the DOWN transition probabilities
    #print deltaT, N, interval, rhoBar, sigmaRhoSquared, tauEff
    UP   = synapticChange.transitionProbability(Npairs,interval,synapticChange.rhoStar,0.,rhoBar,sigmaRhoSquared,tauEff)
    DOWN = synapticChange.transitionProbability(Npairs,interval,synapticChange.rhoStar,1.,rhoBar,sigmaRhoSquared,tauEff)
    
    # mean value of the synaptic strength right at the end of the stimulation protocol
    meanUP   =  rhoBar - rhoBar*exp(-Npairs*interval/tauEff)
    meanDOWN =  rhoBar - (rhoBar - 1.)*exp(-Npairs*interval/tauEff)
    
    # change in synaptic strength after/before
    synChange = synapticChange.changeSynapticStrength(synapticChange.beta,UP,DOWN,synapticChange.b)
    
    return synChange

#######################################################

#try:
#    paramSets = array([int(sys.argv[1])])
#    #oldNewSets = str(sys.argv[2])
#except IndexError:
#    paramSets  = array([0,2,7,9,10,24,27,29,37])
    #oldNewSets = array(['new','new','new','new','new','new','new','new','new'])

#print paramSets
#print oldNewSets
#pdb.set_trace()
#dataSetNumber = 37

Ca0 = 2.
CaTest = 2.

#######################################################
# load experimental data
dataLocation = 'experimental_data/'
stdp1Hz = np.loadtxt(dataLocation+'STDP_1Hz_100pairings.dat')
stdp3Hz = np.loadtxt(dataLocation+'STDP_2.5-3Hz_100pairings.dat')
stdp5Hz = np.loadtxt(dataLocation+'STDP_5Hz_100pairings.dat')
stdp10Hz = np.loadtxt(dataLocation+'STDP_10Hz_100pairings.dat')

# reading the binned data
stdp1binned = np.loadtxt(dataLocation+'STDP_1Hz_100pairings_binned.dat')
stdp3binned = np.loadtxt(dataLocation+'STDP_2.5-3Hz_100pairings_binned.dat')
stdp5binned = np.loadtxt(dataLocation+'STDP_5Hz_100pairings_binned.dat')
stdp10binned = np.loadtxt(dataLocation+'STDP_10Hz_100pairings_binned.dat')

stimFrequencies = [1,3,5,10]

#########################################################
# load fit results
#fitLocation = '/home/mgraupe/theobio/camkII/simplest_model/laurent_venance/parameter_search/output/summary_check_outcome/'
#fit = np.loadtxt(fitLocation+'venance_transition_outcome_ds#0.dat')

amountOfChange = zeros((len(paramSets),3))

dsN = 0

solutions = pickle.load(open('solutions.py'))
print solutions

for dataSetNumber in paramSets:
    #oldNew = oldNewSets[dsN] 
    print dataSetNumber #, oldNew
    #########################################################
    # analytical results for irregular pairs
    
    #if 'old' in oldNew:
    #    analyticalLocation = '/home/mgraupe/theobio/camkII/simplest_model/laurent_venance/irregular_pairs/parameterset#%s/' % dataSetNumber
    #    model1Hz = np.loadtxt(analyticalLocation+'freq1Hz/poisson_transition_prob.dat')
    #    model3Hz = np.loadtxt(analyticalLocation+'freq3Hz/poisson_transition_prob.dat')
    #    model5Hz = np.loadtxt(analyticalLocation+'freq5Hz/poisson_transition_prob.dat')
    #    model10Hz = np.loadtxt(analyticalLocation+'freq10Hz/poisson_transition_prob.dat')

    #if 'new' in oldNew:
    #analyticalLocation = '/home/mgraupe/theobio/camkII/simplest_model/laurent_venance/irregular_pairs/parameterset#%s/' % dataSetNumber
    #modelNew     = np.loadtxt(analyticalLocation+'irregular_spike_pairs.dat')
    #modelNewFreq = np.loadtxt(analyticalLocation+'irregular_spike_pairs_frequency.dat')
    
    # reverse erronous mirroing of the matrix 
    #modelNewFreq[:,0] = modelNewFreq[:,0][::-1]
    #modelNewFreq = modelNewFreq[::-1]
    
    #pdb.set_trace()
    ##########################################################
    synapticChange.choseParameterSet(analyticalLocation+'parameters.par',fromFile=True)

    ##########################################################
    # generate analytical results for regular pairs
    Npairs = 100.
    tat = timeAboveThreshold(synapticChange.thetaD, synapticChange.thetaP, synapticChange.tauCa, synapticChange.Cpre*CaTest/Ca0, synapticChange.Cpost*CaTest/Ca0)
    
    ###########################################################
    # change as function of delta t for four different frequencies
    deltaTstart=-0.3
    deltaTend  =0.3
    steps = 3001.

    deltaT = linspace(deltaTstart,deltaTend,steps)
    synChange = zeros((len(deltaT),5))
    
    frequency = 1.
    interval = 1./frequency
    for n in range(len(deltaT)):
        synChange[n,0] = deltaT[n]
        synChange[n,1] = calculateSnapticeChange(deltaT[n],frequency,interval,tat)
    #
    frequency = 3.
    interval = 1./frequency
    for n in range(len(deltaT)):
        synChange[n,2] = calculateSnapticeChange(deltaT[n],frequency,interval,tat)
    #
    frequency = 5.
    interval = 1./frequency
    for n in range(len(deltaT)):
        synChange[n,3] = calculateSnapticeChange(deltaT[n],frequency,interval,tat)
    #
    frequency = 10.
    interval = 1./frequency
    for n in range(len(deltaT)):
        synChange[n,4] = calculateSnapticeChange(deltaT[n],frequency,interval,tat)
    
    # ###########################################################
    # # change as function of frequency for four different delta t values
    #
    # freqStart = 0.1
    # freqEnd   = 10.
    # freqSteps = 2001
    #
    # frequencies = linspace(freqStart,freqEnd,freqSteps)
    # synChangeFreq = zeros((len(frequencies),5))
    #
    # deltaT = 0.05
    # for n in range(len(frequencies)):
    #     synChangeFreq[n,0] = frequencies[n]
    #     synChangeFreq[n,1] = calculateSnapticeChange(deltaT,frequencies[n],1./frequencies[n],tat)
    # #
    # deltaT = 0.1
    # for n in range(len(frequencies)):
    #     synChangeFreq[n,2] = calculateSnapticeChange(deltaT,frequencies[n],1./frequencies[n],tat)
    # #
    # deltaT = 0.15
    # for n in range(len(frequencies)):
    #     synChangeFreq[n,3] = calculateSnapticeChange(deltaT,frequencies[n],1./frequencies[n],tat)
    # #
    # deltaT = 0.2
    # for n in range(len(frequencies)):
    #     synChangeFreq[n,4] = calculateSnapticeChange(deltaT,frequencies[n],1./frequencies[n],tat)
    
    
    #pdb.set_trace()
    # mask = (synChange[:,0]*1000.>15.) & (synChange[:,0]*1000.<20.)
    # amountOfChange[dsN,0] = synapticChange.tauCa
    # amountOfChange[dsN,1] = mean(synChange[:,1][mask]*100.)
    
    #######################################################
    # plot data
    fig_width = 12 # width in inches
    fig_height = 10  # height in inches
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
            'xtick.major.size' : 4,      # major tick size in points
            #'edgecolor' : None
            #'xtick.major.size' : 2,
            #'ytick.major.size' : 2,
            }
    rcParams.update(params)
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #plt.rc('text', usetex=True)
    #plt.rcParams['mathtext.fontset'] = 'custom'
    #plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    #plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    #plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    #plt.rc('text', usetex=True)
    # set sans-serif font to Arial
    #rcParams['font.sans-serif'] = 'Arial'
    #plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='sans-serif')
    
    # create figure instance
    fig = plt.figure()


    # define sub-panel grid and possibly width and height ratios
    gs = gridspec.GridSpec(2, 2
                        #width_ratios=[1,1.2],
                        #height_ratios=[1,1]
                        )

    # define vertical and horizontal spacing between panels
    gs.update(wspace=0.2,hspace=0.35)

    fig.suptitle(r'STDP data, data-set %s, $C_{\rm pre} = %s$, $C_{\rm post} = %s$, (100 regular- and irregular-spaced spike-pairs)' % (dataSetNumber, synapticChange.Cpre, synapticChange.Cpost) ,fontsize=14)
    # possibly change outer margins of the figure
    plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.08)

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
    #ax0.axvline(x=-50,ls='--',color='turquoise',lw=2)
    #ax0.axvline(x=25,ls='--',color='turquoise',lw=2)
    #ax0.axvline(x=-200,ls='--',color='turquoise',lw=2)
    #ax0.axvline(x=200,ls='--',color='turquoise',lw=2)
    ax0.plot(stdp1Hz[:,0],stdp1Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
    ax0.errorbar(stdp1binned[:,0],stdp1binned[:,1],yerr=stdp1binned[:,2],fmt='o-',markeredgecolor='C0')
    ax0.plot(synChange[:,0]*1000.,synChange[:,1]*100.)
    #if oldNew == 'old':
    #    #print len(-model1Hz[:,3]+2.*synapticChange.D*1000.), -model1Hz[:,3]+2.*synapticChange.D*1000.
    #    ax0.plot(-model1Hz[:,3]+2.*synapticChange.D*1000.,model1Hz[:,18]*100.)
    #    
    #    mask = ((-model1Hz[:,3]+2.*synapticChange.D*1000.)>15.) & ((-model1Hz[:,3]+2.*synapticChange.D*1000.)<20.)
    #    amountOfChange[dsN,2] = mean(model1Hz[mask,18]*100.)
        
    #else:
    #print len(modelNew[:,0]*1000.), modelNew[:,0]*1000.
    ax0.plot(modelNew[:,0]*1000.,modelNew[:,1]*100.)
    #ax0.plot(-50,modelNewFreq[4,5]*100.,'o')

    mask = (modelNew[:,0]*1000.>15.) & (modelNew[:,0]*1000.<20.)
    amountOfChange[dsN,2] = mean(modelNew[mask,1]*100.)
    
    
    # removes upper and right axes 
    # and moves left and bottom axes away
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_position(('outward', 10))
    ax0.spines['left'].set_position(('outward', 10))
    ax0.yaxis.set_ticks_position('left')
    ax0.xaxis.set_ticks_position('bottom')

    ax0.set_xlim(-300,300)
    # legends and labels
    #plt.legend(loc=1,frameon=False)

    #plt.xlabel(r'$\Delta t$ (ms)')
    plt.ylabel('change in synaptic strength')


    ax_inset = plt.axes((0.35, 0.71, 0.15, 0.2))
    #plt.hist(residuals, fc='0.8',ec='w', lw=2)
    ax_inset.axhline(y=100,ls='--',color='0.7',lw=2)
    ax_inset.axvline(x=0,ls='--',color='0.7',lw=2)
    #ax_inset.axvline(x=-50,ls='--',color='turquoise',lw=2)
    #ax_inset.axvline(x=25,ls='--',color='turquoise',lw=2)
    ax_inset.plot(stdp1Hz[:,0],stdp1Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
    ax_inset.errorbar(stdp1binned[:,0],stdp1binned[:,1],yerr=stdp1binned[:,2],fmt='o-',markeredgecolor='C0')
    ax_inset.plot(synChange[:,0]*1000.,synChange[:,1]*100.)
    #if oldNew == 'old':
    #    ax_inset.plot(-model1Hz[:,3]+2.*synapticChange.D*1000.,model1Hz[:,18]*100.)
    #else:
    ax_inset.plot(modelNew[:,0]*1000.,modelNew[:,1]*100.)
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
    #ax1.axvline(x=-100,ls='--',color='turquoise',lw=2)
    #ax1.axvline(x=100,ls='--',color='turquoise',lw=2)
    #ax1.axvline(x=-250,ls='--',color='turquoise',lw=2)
    #ax1.axvline(x=250,ls='--',color='turquoise',lw=2)
    ax1.plot(stdp3Hz[:,0],stdp3Hz[:,1],'o',ms=4,c='0.5',markeredgecolor='0.5')
    ax1.errorbar(stdp3binned[:,0],stdp3binned[:,1],yerr=stdp3binned[:,2],fmt='o-',markeredgecolor='C0',label=r'exp. data: reg. pairs, Ca$_{\rm ext} = 2$ mM')
    ax1.plot(synChange[:,0]*1000.,synChange[:,2]*100.,label=r'model fit: reg. pairs, Ca$_{\rm ext} = 2$ mM')
    #if 'old'==oldNew:
    #ax1.plot(-model3Hz[:,3]+2.*synapticChange.D*1000.,model3Hz[:,18]*100.)
    #else:
    ax1.plot(modelNew[:,0]*1000.,modelNew[:,2]*100.,label=r'model prediction: irr. pairs, Ca$_{\rm ext} = 2$ mM')

    # removes upper and right axes 
    # and moves left and bottom axes away
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.spines['left'].set_position(('outward', 10))
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    ax1.set_xlim(-300,300)
    # legends and labels
    plt.legend(frameon=False,loc=1)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    
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
    ax2.errorbar(stdp5binned[:,0],stdp5binned[:,1],yerr=stdp5binned[:,2],fmt='o-',markeredgecolor='C0')
    ax2.plot(synChange[:,0]*1000.,synChange[:,3]*100.)
    #if 'old'==oldNew:
    #ax2.plot(-model5Hz[:,3]+2.*synapticChange.D*1000.,model5Hz[:,18]*100.)
    #else:
    ax2.plot(modelNew[:,0]*1000.,modelNew[:,3]*100.)

    # removes upper and right axes 
    # and moves left and bottom axes away
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_position(('outward', 10))
    ax2.spines['left'].set_position(('outward', 10))
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')

    ax2.set_xlim(-300,300)
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
    ax3.errorbar(stdp10binned[:,0],stdp10binned[:,1],yerr=stdp10binned[:,2],fmt='o-',markeredgecolor='C0')
    ax3.plot(synChange[:,0]*1000.,synChange[:,4]*100.)
    #if 'old'==oldNew:
    #ax3.plot(-model10Hz[:,3]+2.*synapticChange.D*1000.,model10Hz[:,18]*100.)
    #else:
    ax3.plot(modelNew[:,0]*1000.,modelNew[:,4]*100.)

    # removes upper and right axes 
    # and moves left and bottom axes away
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_position(('outward', 10))
    ax3.spines['left'].set_position(('outward', 10))
    ax3.yaxis.set_ticks_position('left')
    ax3.xaxis.set_ticks_position('bottom')

    ax3.set_xlim(-300,300)
    # legends and labels
    #plt.legend(loc=1,frameon=False)

    plt.xlabel(r'$\Delta t$ (ms)')
    #plt.ylabel('change in synaptic strength')
    
    
       
    ## save figure ############################################################
    fname = os.path.basename(__file__)[:-3] + '_dataSet#' + str(dataSetNumber)

    savefig(fname+'.png')
    savefig(fname+'.pdf')
    clf()
    dsN+=1




#######################################################
# plot data
fig_width = 6 # width in inches
fig_height = 12  # height in inches
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
gs = gridspec.GridSpec(3, 1
                    #width_ratios=[1,1.2],
                    #height_ratios=[1,1]
                    )

# define vertical and horizontal spacing between panels
gs.update(wspace=0.2,hspace=0.3)

fig.suptitle(r'STDP data, Laurent Venance lab, LTD between 15 and 20 ms' ,fontsize=14)
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
ax0.set_title(r'$\tau_{\rm Ca}$ vs. regular pair avg. plasticity between 15 and 20 ms')

# diplay of data
ax0.axhline(y=100,ls='--',color='0.7',lw=2)
#ax0.axvline(x=0,ls='--',color='0.7',lw=2)
ax0.plot(amountOfChange[:,0]*1000.,amountOfChange[:,1],'o')


# removes upper and right axes 
# and moves left and bottom axes away
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_position(('outward', 10))
ax0.spines['left'].set_position(('outward', 10))
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('bottom')
plt.ylabel(r'$w_{\rm regular}$')


# panel 1 #############################################
ax1 = plt.subplot(gs[1])

# title
ax1.set_title(r'$\tau_{\rm Ca}$ vs. irregular pair avg. plasticity between 15 and 20 ms')

# diplay of data
ax1.axhline(y=100,ls='--',color='0.7',lw=2)
ax1.plot(amountOfChange[:,0]*1000.,amountOfChange[:,2],'o')

# removes upper and right axes 
# and moves left and bottom axes away
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_position(('outward', 10))
ax1.spines['left'].set_position(('outward', 10))
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')
plt.ylabel(r'$w_{\rm irregular}$')


# panel 1 #############################################
ax2 = plt.subplot(gs[2])

# title
ax2.set_title(r'$\tau_{\rm Ca}$ vs. difference in avg. plasticity between 15 and 20 ms')

# diplay of data
#ax2.axhline(y=100,ls='--',color='0.7',lw=2)
ax2.plot(amountOfChange[:,0]*1000.,amountOfChange[:,2]-amountOfChange[:,1],'o')

# removes upper and right axes 
# and moves left and bottom axes away
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_position(('outward', 10))
ax2.spines['left'].set_position(('outward', 10))
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')
plt.ylabel(r'$w_{\rm irregular} - w_{\rm regular}$')

plt.xlabel(r'$\tau_{\rm Ca}$ (ms)')
## save figure ############################################################
fname = os.path.basename(__file__)[:-3] 

savefig(fname+'_tauCa.png')
savefig(fname+'_tauCa.pdf')
