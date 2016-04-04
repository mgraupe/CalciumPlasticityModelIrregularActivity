from pylab import *
import os
import pickle

from synUtils import *
import params as par


def readSol(filenN):
    return pickle.load(open(filenN))

##############################################################################
# instance of synaptic Change and figure class  
synU = synUtils(par.thetaD,par.thetaP,par.nonlinear,par.Npresentations,par.w0)

s1 = readSol('solutions.py')

sorted(s1, key=lambda s1: s1[1])

#synU.generateFig()