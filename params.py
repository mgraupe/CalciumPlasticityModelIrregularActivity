from collections import *

# fixed parameter and parameter ranges
thetaP = 1.5 #2.009289 #1.3 #1.3
nonlinear = 1. #2.
w0 = 0.5
modelVersion =  'multiplicative' #''additive' , 'multiplicative'
fitWeights = [2.,1.,1.,1.]

# fitting loops 
Nruns = 100
threshold = 5. # lower threshold to include parameter set

# smoothnessWeight = 1.

# parameterLimits
limits = OrderedDict([
    ('tauCa',[0.02,0.2]),
    ('Cpre',[1.,2.]),
    ('Cpost',[0.1,1.]),
    #('thetaP',[1.2,4.]),
    ('gammaD',[0.1,1000.]),
    ('gammaP',[0.1,1000.]),
    ('sigma',[0.1,100.]),
    ('tau',[1.,50000.]),
    ('D',[-0.01,0.015]),
    ('b',[1.,10.]),
])

