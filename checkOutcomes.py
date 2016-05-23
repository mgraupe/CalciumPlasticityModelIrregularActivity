from pylab import *
import os
import pickle

from synUtils import *
import params as par


def readSol(filenN):
    sol = pickle.load(open(filenN))
    return sorted(sol, key=lambda sol: sol[1])

def sort(sol):
    return sorted(sol, key=lambda sol: sol[1])

def generateAllFigs(sol,fn):
    for i in range(len(sol)):
        print i
        synU.generateFig(sol[i],figName=fn+str(i))



##############################################################################
# instance of synaptic Change and figure class  
synU = synUtils(par.thetaD,par.thetaP,par.nonlinear,par.Npresentations,par.w0)

#s5 = readSol('solutions_stoch5.py')
#s6 = readSol('solutions_stoch6.py')
s0 = readSol('solutions_thetap0.py')
s1 = readSol('solutions_thetap1.py')

#generateAllFigs(s13,'solutionsStoch13_')
#generateAllFigs(s14,'solutionsStoch14_')

#sorted(s1, key=lambda s1: s1[1])

#synU.generateFig()