from collections import *

# fixed parameter and parameter ranges
thetaD = 1.
thetaP = 1.3 #1.3
nonlinear = 1. #2.
Npresentations = 75.
w0 = 0.5

# fitting loops 
Nruns = 10
threshold = 1000. # lower threshold to include parameter set

# parameterLimits
limits = OrderedDict([
    ('tauCa',[0.001,0.1]),
    ('Cpre',[0.01,1.3]),
    ('Cpost',[0.01,1.3]),
    ('gammaD',[0.1,1000.]),
    ('gammaP',[0.1,1000.]),
    ('tau',[1.,50000.]),
    ('D',[0.,0.05]),
])
