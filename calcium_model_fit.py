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
def errFunc(params, stimFrequencies, fitWeights, fitData1Hz,fitData3Hz,fitData5Hz,fitData10Hz):
    # mse for data vs. model
    chi2 = 0.
    for n in range(len(stimFrequencies)):
        exec('dataSet = fitData%sHz' % stimFrequencies[n])
        # print stimFrequencies[n], dataSet

        for i in range(len(dataSet)):
                yModel = synU.calculateChangeInSynapticStrength(stimFrequencies[n],dataSet[i,0]/1000.,params)
                #
                # print yModel, dataSet[i,1]/100.
                chi2+= fitWeights[n]*((dataSet[i,1]/100. - yModel)**2)/((dataSet[i,2]/100.)**2)

    # add smoothness constraint
    #frequencies = linspace(1.,50.,50)
    #synChangePlus  = zeros(len(frequencies))
    #synChangeMinus = zeros(len(frequencies))
    #for n in range(len(frequencies)):
    #    synChangePlus[n]  = synU.calculateChangeInSynapticStrength(frequencies[n],0.01,params)
    #    synChangeMinus[n] = synU.calculateChangeInSynapticStrength(frequencies[n],-0.01,params)
    #chi2 += par.smoothnessWeight*(sum((mean(synChangePlus)-synChangePlus)**2))
    #chi2 += par.smoothnessWeight*(sum((mean(synChangeMinus)-synChangeMinus)**2))
    # penalize if parameters are outside pre-defined limits
    i=0
    for k in par.limits:
        #print k, params[i]
        if (params[i] < par.limits[k][0]) or (params[i] > par.limits[k][1]):
            #print 'penalized'
            chi2+=100.
        i+=1
    # impose Cpost>Cpre
    #if params[1] > params[2]:
    #    chi2+=100.
    
    return chi2


############################################################################
def initialGuess(lim):
    nParams = len(lim)
    params = zeros(nParams)
    random.seed(int64((time.time()-base)*100))
    randTemp = rand(nParams)
    n = 0
    for k in lim:
        params[n] = lim[k][0] + randTemp[n]*(lim[k][1]-lim[k][0])
        #exec('k = %s' % initP)
        n+=1
    return params


#os.nice(19)
base = 1531818909

##############################################################################
# instance of synaptic Change and figure class
synU = synUtils(par.thetaD,par.thetaP,par.nonlinear,par.Npresentations,par.w0,dataSet='venance')

#############################################################################
#pdb.set_trace()
solutions = []
for n in range(par.Nruns):
    #Initial guess of parameters
    #params0  = [0.0667179, 1.45248, 0.405039, 2.0, 15.973, 16.3457, -0.00156591] #
    params0 = initialGuess(par.limits)
    #pdb.set_trace()
    # Apply downhill Simplex algorithm.
    p1 = simplex(errFunc, params0, args=(synU.stimFrequencies,synU.fitWeights,synU.rawData1Hz,synU.rawData3Hz,synU.rawData5Hz,synU.rawData10Hz), full_output=1, disp=True,maxiter=1E4, maxfun=1E4)
    #print p1
    if p1[1] < par.threshold: 
        solutions.append(p1)

if solutions:
    sorted(solutions, key=lambda solutions: solutions[1])
    pickle.dump(solutions,open('solutions.py','w'))




