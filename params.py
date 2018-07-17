from collections import *

# fixed parameter and parameter ranges
thetaD = 1.
thetaP = 1.3 #1.3
nonlinear = 1. #2.
Npresentations = 100.
w0 = 0.5

# fitting loops 
Nruns = 100000
threshold = 100. # lower threshold to include parameter set

# smoothnessWeight = 1.

# parameterLimits
limits = OrderedDict([
    ('tauCa',[0.001,0.1]),
    ('Cpre',[0.01,0.8]),
    ('Cpost',[0.8,1.3]),
    ('gammaD',[0.1,1000.]),
    ('gammaP',[0.1,1000.]),
    ('tau',[1.,50000.]),
    ('D',[-0.005,0.005]),
])

