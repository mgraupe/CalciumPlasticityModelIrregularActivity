from pylab import *
import os
import pickle

from synUtils import *
import params as par
import parameter_fit_solutions as pfs


def readSol(filenN):
    sol = pickle.load(open(filenN))
    return sorted(sol, key=lambda sol: sol[1])

def sort(sol):
    return sorted(sol, key=lambda sol: sol[1])

def generateAllFigs(sol,fn):
    for i in range(len(sol)):
        print(i)
        synU.generateFig(sol[i],figName=fn+str(i))




##############################################################################
# instance of synaptic Change and figure class  
synU = synUtils(dataSet='venance',modelV=par.modelVersion)

#s0 = readSol('solutions.p')
#s6 = readSol('solutions_stoch6.py')
#s0 = readSol('solutions_thetap0.py')
#s1 = readSol('solutions_thetap1.py')

#synU.generateVenFig(pfs.VenanceSmult5,modelV='m')
#synU.generateVenFig('VenanceSmult9',modelV='m')
#synU.generateVenFig('linearCaModel',modelV='m')
synU.generateVenFig('VenancesBin0',modelV='bin')
#synU.generateVenFig('VenanceSmult7',modelV='m')
#generateAllFigs(s13,'solutionsStoch13_')
#generateAllFigs(s14,'solutionsStoch14_')


#synU.generateFig()