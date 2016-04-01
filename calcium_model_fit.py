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

def errFunc(params, xData, yData, errData):
    # compute chi-square
    chi2 = 0.
    for n in range(len(xData)):
        yModel = synU.calculateChangeInSynapticStrength(xData[n,0],xData[n,1],params)
        #
        chi2+= (yData[n] - yModel)*(yData[n] - yModel) #/(errData[n]*errData[n])
    i=0
    for k in par.limits:
        if (params[i] < par.limits[k][0]) or (params[i] > par.limits[k][1]):
           chi2+=100.
        i+=1
    return chi2

def initialGuess(limits):
    nParams = len(limits)
    params = zeros(nParams)
    randTemp = rand(nParams)
    n = 0
    for k in limits:
        params[n] = limits[k][0] + randTemp[n]*(limits[k][1]-limits[k][0])
        #exec('k = %s' % initP)
        n+=1
    return params


##############################################################################
# instance of synaptic Change and figure class
synU = synUtils(par.thetaD,par.thetaP,par.nonlinear,par.Npresentations,par.w0)

#############################################################################

# read in experimental data
dataDir = 'experimental_data/'

jesperReg = loadtxt(dataDir+'sjoestroem_regular_all.dat')
jesperStoch = loadtxt(dataDir+'sjoestroem_stochastic.dat')
#read_datafile("../experimental_data/sjoestroem_regular_all.dat",jesper_regular,0);
#	// experimental data - random delta-t ( scan_par[x][2] is now the error)
#  	//read_datafile("../experimental_data/sjoestroem_stochastic.dat",jesper_stochastic,1);
xData = jesperReg[:,[0,1]]
xData[:,1] = xData[:,1]/1000. # everything in sec
yData = jesperReg[:,2]+1. # Sjoestroem's data is normalized to 0
sigmaData = jesperReg[:,3]

solutions = []
for n in range(par.Nruns):
    #Initial guess of parameters
    x0  = initialGuess(par.limits)

    # Apply downhill Simplex algorithm.
    p1 = simplex(errFunc, x0, args=(xData, yData, sigmaData), full_output=1, disp=True,maxiter=1E4, maxfun=1E4)
    
    if p1[1] < 0.3: 
        solutions.append(p1)
    #errFunc(p1[0],xData,yData,sigmaData)

pickle.dump(solutions,open('solutions.py','w'))




