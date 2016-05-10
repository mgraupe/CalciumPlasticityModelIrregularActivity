from pylab import *
from scipy.optimize import fmin as simplex
import os
import pickle

from synUtils import *
import params as par

# "fmin" is not a sensible name for an optimisation package.
# Rename fmin to "simplex"

# Define the objective function to be minimised by Simplex.
# params ... array holding the values of the fit parameters.
# X      ... array holding x-positions of observed data.
# Y      ... array holding y-values of observed data.
# Err    ... array holding errors of observed data.

###########################################################################
# the calcuation of mse is defined here 
def errFunc(params, xDataReg, yDataReg, errData, xDataStoch, yDataStoch, errStochData):
    # mse for regular Sjoestroem data
    chi2 = 0.
    for n in range(len(xDataReg)):
        yModel = synU.calculateChangeInSynapticStrength(xDataReg[n,0],xDataReg[n,1],params)
        #
        chi2+= (yDataReg[n] - yModel)*(yDataReg[n] - yModel)/(errData[n]*errData[n])
    # mse for stochastic Sjoestroem data
    for n in range(len(xDataStoch)):
        yModel = synU.calculateChangeInSynapticStrengthStochastic(xDataStoch[n],params,[-0.015,0.015])
        #
        #print n, xDataStoch[n]
        chi2+= (yDataStoch[n] - yModel)*(yDataStoch[n] - yModel)/(errStochData[n]*errStochData[n])
    # add smoothness constraint
    frequencies = linspace(5.,50.,10.)
    synChangePlus  = zeros(len(frequencies))
    synChangeMinus = zeros(len(frequencies))
    for n in range(len(frequencies)):
        synChangePlus[n]  = synU.calculateChangeInSynapticStrength(frequencies[n],0.01,params)
        synChangeMinus[n] = synU.calculateChangeInSynapticStrength(frequencies[n],-0.01,params)
    chi2 += par.smoothnessWeight*(sum((mean(synChangePlus)-synChangePlus)**2))
    chi2 += par.smoothnessWeight*(sum((mean(synChangeMinus)-synChangeMinus)**2))
    # penalize if parameters are outside pre-defined limits
    i=0
    for k in par.limits:
        if (params[i] < par.limits[k][0]) or (params[i] > par.limits[k][1]):
           chi2+=100.
        i+=1
    return chi2


############################################################################
def initialGuess(lim):
    nParams = len(lim)
    params = zeros(nParams)
    randTemp = rand(nParams)
    n = 0
    for k in lim:
        params[n] = lim[k][0] + randTemp[n]*(lim[k][1]-lim[k][0])
        #exec('k = %s' % initP)
        n+=1
    return params


##############################################################################
# instance of synaptic Change and figure class
synU = synUtils(par.thetaD,par.thetaP,par.nonlinear,par.Npresentations,par.w0)

#############################################################################

solutions = []
for n in range(par.Nruns):
    #Initial guess of parameters
    params0  = initialGuess(par.limits)

    # Apply downhill Simplex algorithm.
    p1 = simplex(errFunc, params0, args=(synU.xDataReg, synU.yDataReg, synU.sigmaDataReg, synU.xDataStoch, synU.yDataStoch, synU.sigmaDataStoch), full_output=1, disp=True,maxiter=1E4, maxfun=1E4)
    
    if p1[1] < par.threshold: 
        solutions.append(p1)

if solutions:
    pickle.dump(solutions,open('solutions.py','w'))




