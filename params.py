from collections import *

# fixed parameter and parameter ranges
thetaD = 1.
thetaP = 1.3 #1.3
nonlinear = 1. #2.
Npresentations = 100.
w0 = 0.5
modelVersion =  'multiplicative' #''additive'
fitWeights = [5.,1.,1.,1.]

# fitting loops 
Nruns = 10
threshold = 300. # lower threshold to include parameter set

# smoothnessWeight = 1.

# parameterLimits
limits = OrderedDict([
    ('tauCa',[0.01,0.1]),
    ('Cpre',[1.2,2.]),
    ('Cpost',[0.1,1.]),
    ('gammaD',[0.1,1000.]),
    ('gammaP',[0.1,1000.]),
    ('tau',[100.,50000.]),
    ('D',[-0.005,0.015]),
])

