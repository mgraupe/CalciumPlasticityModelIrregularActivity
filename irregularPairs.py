'''
        Script to calculte the change in synaptic strength for irregular spike-pairs. 
        
        Model defined in :
        Graupner M and Brunel N (2012). 
        Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location. 
        PNAS 109 (10): 3991-3996.
        
'''

from scipy import *
from numpy import *
from pylab import *
import os
import sys
import time
import math
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import multiprocessing
import pdb
#import commands

#import synapticChange
import params as par
#from synUtils import *
import bestSolutions as bS
from timeAboveThreshold import timeAboveThreshold
from synapticChange import synapticChange

##########################################################
def runIrregularPairSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    #if synChange.Cpre > synChange.Cpost:
    #    (alphaD, alphaP) = tat.irregularSpikePairsEventBased(dT + synChange.D, preRate, postRate, p)
    #else:
    (alphaD, alphaP) = tat.irregularSpikePairsEventBased(dT - synChange.D, preRate, postRate, p)
    #synChange.changeInSynapticStrength(N_pres/preRate, rho0, alphaD, alphaP)
    synChange.changeInSynapticStrength(N_pres/preRate, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
def runRegularPairSimulations(args):
    dT = args[0]
    preRate = args[1]
    postRate = args[2]
    p = args[3]

    (alphaD, alphaP) = tat.spikePairFrequencyNonlinear(dT - synChange.D, preRate)
    #print dT, preRate, alphaD, alphaP
    synChange.changeInSynapticStrength(N_pres/preRate, rho0, alphaD, alphaP)

    return synChange.mean


##########################################################
##########################################################
# output directory
outputDir = 'FigsSimResults/'

# numerical integration step width
#deltaCa = 0.005  # 0.005 #0.01 #  0.0001
#T_total = 100.  # total time of stimulation in sec
N_pres = 100.
rho0 = 0.5
nl = 1.  # nonlinearity factor

###########################################################
# initiate synaptic change class and chose parameter set from file

params = 'VenanceSmult6'
synChange = synapticChange(params, fromFile=True, nonlinear=nl)
# initiate class to calculate fraction of time above threshold
print 'Parameters :',synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP, nl
tat = timeAboveThreshold(synChange.tauCa, synChange.Cpre, synChange.Cpost, synChange.thetaD, synChange.thetaP, nonlinear=nl)

pool = multiprocessing.Pool()

##################################################################################################
# synaptic change vs Delta T for irregular Pairs
#################################################################################################

print 'irregular pairs : synaptic change vs Delta T for four frequencies and one p\'s'

# Parameter of the stimulation protocol
frequencies =  array([1.,3.,5.,10.])  # frequency of spike-pair presentations in pairs/sec
DeltaTstart = -0.3  # start time difference between pre- and post-spike, in sec
DeltaTend = 0.3  # end time difference between pre- and post-spike, in sec
DeltaTsteps = 301.  # steps between start and end value
ppp = 1.

nCases = len(frequencies)

###########################################################
# initialize arrays
deltaT = linspace(DeltaTstart, DeltaTend, DeltaTsteps)
resultsIrr = zeros(len(frequencies) * 3 + 2)

###########################################################
# simulation loop over range of deltaT values
for i in range(len(deltaT)):
    #
    print 'deltaT : ', deltaT[i]

    args = column_stack((ones(nCases) * deltaT[i], frequencies, frequencies, ones(nCases) * ppp))

    rrr = pool.map(runIrregularPairSimulations, args)
    # for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    # pdb.set_trace()
    res1 = hstack((deltaT[i], frequencies, frequencies, ppp, rrr))
    resultsIrr = vstack((resultsIrr, res1))

resultsIrr = resultsIrr[1:]

# the array has to be flipped if Cpre>Cpost
#if synChange.Cpre > synChange.Cpost:
#    # invert time axis
#    resultsIrr[:, 0] = resultsIrr[:, 0][::-1]
#    # reorder entire array to end up with increasing time
#    resultsIrr = resultsIrr[::-1]

if not os.path.exists(outputDir):
    os.makedirs(outputDir)

np.save(outputDir + 'irregularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % params, resultsIrr)
np.savetxt(outputDir + 'irregularSpikePairs_vs_deltaT_differentFreqs_%s.dat' % params, resultsIrr)

##########################################################
# synaptic change vs Delta T for regular Pairs
##########################################################
print 'regular pairs : synaptic change vs Delta T for six frequencies and one p'

resultsReg = zeros(len(frequencies)*3+2)

# simulation loop over range of deltaT values
for i in range(len(deltaT)):
    #
    print 'deltaT : ', deltaT[i]

    args = column_stack((ones(nCases) * deltaT[i], frequencies, frequencies, ones(nCases) * ppp))

    rrr = pool.map(runRegularPairSimulations, args)
    # for n in range(len(frequencies)):
    #    (synC[i,n],meanU[i,n],meanD[i,n],tD[i,n],tP[i,n]) = rrr[n]
    # pdb.set_trace()
    res1 = hstack((deltaT[i], frequencies, frequencies, ppp, rrr))
    resultsReg = vstack((resultsReg, res1))

resultsReg = resultsReg[1:]

np.save(outputDir + 'regularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % params, resultsReg)
np.savetxt(outputDir + 'regularSpikePairs_vs_deltaT_differentFreqs_%s.dat' % params, resultsReg)







