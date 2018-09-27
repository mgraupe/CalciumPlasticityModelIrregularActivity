import numpy as np
from scipy.optimize import fmin as simplex
import pickle
import time
import multiprocessing
import pdb

import params as par

from timeAboveThreshold import timeAboveThreshold
from synapticChange import synapticChange


# "fmin" is not a sensible name for an optimisation package.
# Rename fmin to "simplex"

# Define the objective function to be minimised by Simplex.
# params ... array holding the values of the fit parameters.
# X      ... array holding x-positions of observed data.
# Y      ... array holding y-values of observed data.
# Err    ... array holding errors of observed data.


###########################################################################
def fitData(params):
    print params
    p1 = simplex(errFunc, params, args=(sChange.stimFrequencies,par.fitWeights,sChange.rawData1Hz,sChange.rawData3Hz,sChange.rawData5Hz,sChange.rawData10Hz,sChange.rawIrregularData1Hz,sChange.rawIrregularData3Hz), full_output=1, disp=True,maxiter=1E4, maxfun=1E4)
    return p1

###########################################################################
# the calcuation of mse is defined here 
def errFunc(params, stimFrequencies, fitWeights, fitData1Hz,fitData3Hz,fitData5Hz,fitData10Hz,fitIrrData1Hz,fitIrrData3Hz):
    # mse for data vs. model
    chi2 = 0.
    sChange.choseParameterSet('Venance',source='fromDictionary',params=params)
    # initiate class to calculate fraction of time above threshold
    #print sChange.thetaP
    tat = timeAboveThreshold(sChange.tauCa, sChange.Cpre, sChange.Cpost, sChange.thetaD, sChange.thetaP, nonlinear=par.nonlinear)
    #print 'Parameters :', sChange.tauCa, sChange.Cpre, sChange.Cpost, sChange.thetaD, sChange.thetaP, sChange.D

    for n in range(len(stimFrequencies)):
        exec('dataSet = fitData%sHz' % stimFrequencies[n])
        # print stimFrequencies[n], dataSet
        if fitWeights[n]:
            for i in range(len(dataSet)):
                (alphaD, alphaP) = tat.spikePairFrequencyNonlinear(dataSet[i,0]/1000. - sChange.D, float(stimFrequencies[n]))
                # print dT, preRate, alphaD, alphaP
                sChange.changeInSynapticStrength(sChange.Npresentations/float(stimFrequencies[n]), par.w0, alphaD, alphaP)
                #print float(stimFrequencies[n]), sChange.Npresentations/float(stimFrequencies[n])
                #yModel = synU.calculateChangeInSynapticStrength(float(stimFrequencies[n]),dataSet[i,0]/1000.,params)
                #
                #if stimFrequencies[n] == 1:
                #    chi2 += fitWeights[n] * ((dataSet[i, 1] / 100. - yModel) ** 2) * (dataSet[i, 2])
                #    #print dataSet[i,0]/1000., yModel, dataSet[i,1]/100., dataSet[i,2]
                #else:
                #print dataSet[i,1]/100., sChange.mean/par.w0
                chi2+= fitWeights[n]*((dataSet[i,1]/100. - sChange.mean/par.w0)**2) #*(dataSet[i,2])
                #print sChange.meanAdd
                #chi3+= ((dataSet[i,1]/100. - 1.)**2)/((dataSet[i,2]/100.)**2)
        if stimFrequencies[n] in [1,3]:
            exec('dataSet = fitIrrData%sHz' % stimFrequencies[n])
            #print 'irregular', stimFrequencies[n], dataSet
            for i in range(len(dataSet)): #irregularSpikePairsEventBased(self, deltaT, preRate, postRate, ppp)
                (alphaD, alphaP) = tat.irregularSpikePairsEventBased(dataSet[i,0] - sChange.D, float(stimFrequencies[n]),float(stimFrequencies[n]),1.)
                sChange.changeInSynapticStrength(sChange.Npresentations/float(stimFrequencies[n]), par.w0, alphaD, alphaP)
                chi2+= ((dataSet[i,2] - sChange.mean/par.w0)**2) #*(dataSet[i,2])

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
            #chi3+=100.
        i+=1
    # impose Cpost>Cpre
    #if params[1] > params[2]:
    #    chi2+=100.
    #print chi3
    return chi2


############################################################################
def initialGuess(lim):
    nParams = len(lim)
    params = np.zeros(nParams)
    #np.random.seed(np.int64((time.time()-base)*100))
    randTemp = np.random.rand(nParams)
    n = 0
    for k in lim:
        params[n] = lim[k][0] + randTemp[n]*(lim[k][1]-lim[k][0])
        #exec('k = %s' % initP)
        n+=1
    return params


#os.nice(19)
base = 1537192766

##############################################################################
# instance of synaptic Change and figure class
sChange = synapticChange('Venance') #, source='fromDictionary', nonlinear=par.nonlinear,parameter=params)

#synU = synUtils(par.thetaD,par.thetaP,par.nonlinear,par.Npresentations,par.w0,dataSet='venance',modelV=par.modelVersion)

nProcessors = 5
pool = multiprocessing.Pool(nProcessors)

#############################################################################
#pdb.set_trace()
solutions = []
for n in range(par.Nruns):
    #Initial guess of parameters
    #params0  = [0.0667179, 1.45248, 0.405039, 2.0, 15.973, 16.3457, -0.00156591] #
    params0 = []
    for i in range(nProcessors):
        params0.append(initialGuess(par.limits))
    #pdb.set_trace()
    rrr = pool.map(fitData,params0)

    #pdb.set_trace()
    # Apply downhill Simplex algorithm.
    #print 'Number of free parameters : ', len(params0)
    #p1 = simplex(errFunc, params0, args=(sChange.stimFrequencies,par.fitWeights,sChange.rawData1Hz,sChange.rawData3Hz,sChange.rawData5Hz,sChange.rawData10Hz,sChange.rawIrregularData1Hz,sChange.rawIrregularData3Hz), full_output=1, disp=True,maxiter=1E4, maxfun=1E4)
    #print p1
    for i in range(len(rrr)):
        print rrr[i]
        if rrr[i][1] < par.threshold:
            solutions.append(rrr[i])

    if solutions:
        solutions = sorted(solutions, key=lambda solutions: solutions[1])
        pickle.dump(solutions,open('solutions.py','w'))


print solutions[0]