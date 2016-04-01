from pylab import *
from scipy.optimize import fmin as simplex
from timeAboveThreshold import *
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import *
import os

# "fmin" is not a sensible name for an optimisation package.
# Rename fmin to "simplex"

# Define the objective function to be minimised by Simplex.
# params ... array holding the values of the fit parameters.
# X      ... array holding x-positions of observed data.
# Y      ... array holding y-values of observed data.
# Err    ... array holding errors of observed data.


def calculateChangeInSynapticStrength(frequency,deltaT,params):
    #####
    n=0
    #for k in limits: #range(len(params)):
    #    exec('k = %s' % params[n])
    #    n+=1 
    #pdb.set_trace()
    tauCa = params[0]
    Cpre = params[1]
    Cpost = params[2]
    #thetaD = params[3]
    #thetaP = params[4]
    gammaD = params[3]
    gammaP = params[4]
    tau = params[5]
    D = params[6]
    
    interval    = 1./frequency
    ####
    tat = timeAboveThreshold(thetaD, thetaP, tauCa, Cpre, Cpost, nonlinear)
    (timeD,timeP) = tat.spikePairFrequencyNonlinear(deltaT-D,frequency)
    # average potentiation and depression rates
    GammaP = gammaP*timeP/interval
    GammaD = gammaD*timeD/interval
    # rhoBar: average value of rho in the limit of a very long protocol equivalent to the minimum of the quadratic potentia
    try :
        rhoBar = GammaP/(GammaP + GammaD)
    except RuntimeWarning:
        print GammaP, GammaD
    # tauEff : characteristic time scale of the temporal evolution of the pdf of rho
    tauEff = tau/(GammaP + GammaD)
    #
    # mean value of the synaptic strength right at the end of the stimulation protocol
    mean   =  rhoBar - (rhoBar- 0.5)*exp(-Npresentations*interval/tauEff)
    # change in synaptic strength after/before
    return mean/w0

def errFunc(params, xData, yData, errData):
    # compute chi-square
    chi2 = 0.
    for n in range(len(xData)):
        yModel = calculateChangeInSynapticStrength(xData[n,0],xData[n,1],params)
        #
        chi2+= (yData[n] - yModel)*(yData[n] - yModel) #/(errData[n]*errData[n])
    i=0
    for k in limits:
        if (params[i] < limits[k][0]) or (params[i] > limits[k][1]):
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


def generateFigure(paraOpt):
    
    ####################################
    # calculate solution
    freq = linspace(0.1,50.,500)
    synChange = zeros((len(freq),2))

    deltaTs = linspace(-0.05,0.05,1001)
    synChange2 = zeros(len(deltaTs))

    for i in range(len(freq)):
        #calculateChangeInSynapticStrength(frequency,deltaT,params):
        synChange[i,0] = calculateChangeInSynapticStrength(freq[i],0.01,paraOpt[0])
        synChange[i,1] = calculateChangeInSynapticStrength(freq[i],-0.01,paraOpt[0])
        
    for i in range(len(deltaTs)):
        synChange2[i] = calculateChangeInSynapticStrength(0.1,deltaTs[i],paraOpt[0])

    ###############################################################
    # produce figure


    # set plot attributes

    fig_width = 5 # width in inches
    fig_height = 8  # height in inches
    fig_size =  [fig_width,fig_height]
    params = {'axes.labelsize': 14,
            'axes.titlesize': 13,
            'font.size': 11,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'figure.figsize': fig_size,
            'savefig.dpi' : 600,
            'axes.linewidth' : 1.3,
            'ytick.major.size' : 4,      # major tick size in points
            'xtick.major.size' : 4      # major tick size in points
            #'edgecolor' : None
            #'xtick.major.size' : 2,
            #'ytick.major.size' : 2,
            }
    rcParams.update(params)

    # set sans-serif font to Arial
    rcParams['font.sans-serif'] = 'Arial'

    # create figure instance
    fig = plt.figure()


    # define sub-panel grid and possibly width and height ratios
    gs = gridspec.GridSpec(2, 1,
                        #width_ratios=[1,1.2],
                        #height_ratios=[1,1]
                        )

    # define vertical and horizontal spacing between panels
    gs.update(wspace=0.3,hspace=0.4)

    # possibly change outer margins of the figure
    plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.18)

    # first sub-plot #######################################################
    ax0 = plt.subplot(gs[0])

    # title
    ax0.set_title('regular Sjoestroem , chi2='+str(paraOpt[1]))

    # diplay of data
    ax0.axhline(y=1.,c='0.7')
    ax0.plot(xData[:,0][::2],yData[::2],'s',color='red',clip_on=False)
    ax0.plot(xData[:,0][1::2],yData[1::2],'o',color='blue',clip_on=False)
    ax0.plot(freq,synChange[:,0],color='red')
    ax0.plot(freq,synChange[:,1],color='blue')

    # removes upper and right axes 
    # and moves left and bottom axes away
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_position(('outward', 10))
    ax0.spines['left'].set_position(('outward', 10))
    ax0.yaxis.set_ticks_position('left')
    ax0.xaxis.set_ticks_position('bottom')

    # legends and labels
    plt.legend(loc=1,frameon=False)

    plt.xlabel('frequency (Hz)')
    plt.ylabel('change in synaptic strength')

    # first sub-plot #######################################################
    ax1 = plt.subplot(gs[1])

    # title
    #ax1.set_title('regular Sjoestroem')

    ax1.axhline(y=1.,c='0.7')
    ax1.plot(deltaTs*1000.,synChange2)

    # removes upper and right axes 
    # and moves left and bottom axes away
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.spines['left'].set_position(('outward', 10))
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    # legends and labels
    plt.legend(loc=1,frameon=False)

    plt.xlabel(r'$\Delta t$ (ms)')
    plt.ylabel('change in synaptic strength')

    # change legend text size 
    #leg = plt.gca().get_legend()
    #ltext  = leg.get_texts()
    #plt.setp(ltext, fontsize=11)

    ## save figure ############################################################
    fname = os.path.basename(__file__)

    savefig(fname[:-3]+'.png')
    savefig(fname[:-3]+'.pdf')


##############################################################################
# parameter and parameter ranges
thetaD = 1.
thetaP = 1.3
nonlinear = 2. #2.
Npresentations = 75.
w0 = 0.5
#
Nruns = 1000
# parameterLimits
limits = OrderedDict([
    ('tauCa',[0.001,0.1]),
    ('Cpre',[0.1,1.]),
    ('Cpost',[0.1,1.3]),
    ('gammaD',[0.1,1000.]),
    ('gammaP',[0.1,1000.]),
    ('tau',[1.,50000.]),
    ('D',[0.,0.05]),
])

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
for n in range(Nruns):
    #Initial guess of parameters
    x0  = initialGuess(limits)

    # Apply downhill Simplex algorithm.
    p1 = simplex(errFunc, x0, args=(xData, yData, sigmaData), full_output=1, disp=True,maxiter=1E4, maxfun=1E4)
    
    if p1[1] < 0.3: 
        solutions.append(p1)
    #errFunc(p1[0],xData,yData,sigmaData)

pickle.dump(solutions,open('solutions.py','w'))

#generateFigure(p1[0])



