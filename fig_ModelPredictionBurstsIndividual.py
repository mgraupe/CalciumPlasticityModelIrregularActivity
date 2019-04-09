from pylab import *
import os
import pickle

from synUtils import *
import params as par
import parameter_fit_solutions as pfs


################################################################################################
def generateFitRegDataFig(paraName, dataDir, figDir, modelV=None):

    w0 = 0.5
    exec('paraOpt = pfs.%s' % paraName)
    # synapticChange.choseParameterSet(analyticalLocation + 'parameters.par', fromFile=True)
    stdp1Hz = np.loadtxt(dataDir + 'STDP_1Hz_100pairings.dat')
    stdp3Hz = np.loadtxt(dataDir + 'STDP_2.5-3Hz_100pairings.dat')
    stdp5Hz = np.loadtxt(dataDir + 'STDP_5Hz_100pairings.dat')
    stdp10Hz = np.loadtxt(dataDir + 'STDP_10Hz_100pairings.dat')

    allData = pickle.load(open(dataDir + "allIrregularData.p", "rb")) #, encoding='latin1')

    # irrSTDP1Hz = np.load(self.dataDir + '18-10-05_experimentOverview_1Hz.npy')
    # irrSTDP3Hz = np.load(self.dataDir + '18-10-05_experimentOverview_3Hz.npy')

    ##########################################################
    # generate analytical results for regular pairs
    irrDatDir = 'FigsSimResults/'
    if os.path.isfile(irrDatDir + 'irregularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName):
        irrData = True
        modelNew = np.load(irrDatDir + 'irregularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName)
        modelBurstsNew = np.load(irrDatDir + 'irregularBurstSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName)
        modelIndNew = np.load(irrDatDir + 'irregularIndividualSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName)
        modelNewReg = np.load(irrDatDir + 'regularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName)
    else:
        irrData = False
    # Npairs = 100.
    # tat = timeAboveThreshold(synapticChange.thetaD, synapticChange.thetaP, synapticChange.tauCa, synapticChange.Cpre * CaTest / Ca0, synapticChange.Cpost * CaTest / Ca0)

    ###########################################################
    # change as function of delta t for four different frequencies
    deltaTstart = -0.3
    deltaTend = 0.3
    steps = 3001.

    # stimFreq = [1,3,5,10]

    if modelV is not None:
        if modelV == 'a':
            modelVersion = 'additive'
        elif modelV == 'm':
            modelVersion = 'multiplicative'
        elif modelV == 'bin':
            modelVersion = 'binary'
        else:
            print
            'problem in model choice'
            sys.exit(1)

    sChange = synapticChange('Venance')  # ,threshold=par.thetaP)
    sChange.choseParameterSet(paraName, source='fromFile')
    # initiate class to calculate fraction of time above threshold
    tat = timeAboveThreshold(sChange.tauCa, sChange.Cpre, sChange.Cpost, sChange.thetaD, sChange.thetaP)
    print
    'Parameters :', sChange.tauCa, sChange.Cpre, sChange.Cpost, sChange.thetaD, sChange.thetaP

    deltaT = linspace(deltaTstart, deltaTend, steps)
    synChange = zeros((len(deltaT), len(sChange.stimFrequencies) + 1))
    for n in range(len(sChange.stimFrequencies)):
        # frequency = stimFreq[n]
        interval = 1. / sChange.stimFrequencies[n]
        for i in range(len(deltaT)):
            (alphaD, alphaP) = tat.spikePairFrequencyNonlinear(deltaT[i] - sChange.D, sChange.stimFrequencies[n])
            # print dT, preRate, alphaD, alphaP
            sChange.changeInSynapticStrength(sChange.Npresentations / sChange.stimFrequencies[n], w0, alphaD, alphaP)

            if n == 0:
                synChange[i, 0] = deltaT[i]
            #  calculateChangeInSynapticStrength(self, frequency,deltaT,params):
            # sol = [0.0667179, 1.45248, 0.405039, 2.0, 15.973, 16.3457, -0.00156591]
            if modelVersion == 'additive':
                synChange[i, n + 1] = sChange.meanAdd / w0
            elif modelVersion == 'multiplicative':
                synChange[i, n + 1] = sChange.mean / w0
            elif modelVersion == 'binary':
                synChange[i, n + 1] = sChange.synChange  # print self.stimFrequencies[n], deltaT[i], synChange[i, n+1]
    ##################################################################

    stdp1binned = np.copy(sChange.rawData1Hz)
    stdp3binned = np.copy(sChange.rawData3Hz)
    stdp5binned = np.copy(sChange.rawData5Hz)
    stdp10binned = np.copy(sChange.rawData10Hz)

    stdpIrr1HzBinned = np.copy(sChange.rawIrregularData1Hz)
    stdpIrr3HzBinned = np.copy(sChange.rawIrregularData3Hz)
    # pdb.set_trace()
    #######################################################
    # plot data
    fig_width = 9  # width in inches
    fig_height = 4  # height in inches
    fig_size = [fig_width, fig_height]
    params = {'axes.labelsize': 14, 'axes.titlesize': 13, 'font.size': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'figure.figsize': fig_size,  # 'savefig.dpi': 600,
              'axes.linewidth': 1.3, 'ytick.major.size': 4,  # major tick size in points
              'xtick.major.size': 4,  # major tick size in points
              # 'edgecolor' : None
              # 'xtick.major.size' : 2,
              # 'ytick.major.size' : 2,
              }
    rcParams.update(params)
    # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # plt.rc('text', usetex=True)
    # plt.rcParams['mathtext.fontset'] = 'custom'
    # plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    # plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    # plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

    # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})
    # plt.rc('text', usetex=True)
    # set sans-serif font to Arial
    # rcParams['font.sans-serif'] = 'Arial'
    # plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='sans-serif')

    # create figure instance
    fig = plt.figure()

    # define sub-panel grid and possibly width and height ratios
    gs = gridspec.GridSpec(1, 2  # width_ratios=[1,1.2],
                           # height_ratios=[1,1]
                           )

    # define vertical and horizontal spacing between panels
    gs.update(wspace=0.3, hspace=0.25)

    #fig.suptitle(r'STDP data, %s, $C_{\rm pre} = %s$, $C_{\rm post} = %s$, (RMS = %s)' % (modelVersion, np.round(sChange.Cpre, 4), np.round(sChange.Cpost, 4), np.round(paraOpt[1], 4)),
    #    fontsize=14)
    # possibly change outer margins of the figure
    plt.subplots_adjust(left=0.1, right=0.96, top=0.92, bottom=0.16)

    # sub-panel enumerations
    plt.figtext(0.0, 0.93, 'A',clip_on=False,color='black', weight='bold',size=22)
    plt.figtext(0.51, 0.93, 'B',clip_on=False,color='black', weight='bold',size=22)
    #plt.figtext(0.0, 0.47, 'C',clip_on=False,color='black', weight='bold',size=22)
    #plt.figtext(0.51, 0.47, 'D',clip_on=False,color='black', weight='bold',size=22)

    # panel 0 #######################################################
    ax0 = plt.subplot(gs[0])

    # title
    ax0.set_title('1 Hz')

    # diplay of data
    ax0.axhline(y=100, ls='--', color='0.7', lw=2)
    ax0.axvline(x=0, ls='--', color='0.7', lw=2)
    # ax0.axvline(x=-50,ls='--',color='turquoise',lw=2)
    # ax0.axvline(x=25,ls='--',color='turquoise',lw=2)
    # ax0.axvline(x=-200,ls='--',color='turquoise',lw=2)
    # ax0.axvline(x=200,ls='--',color='turquoise',lw=2)
    #ax0.plot(stdp1Hz[:, 0], stdp1Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
    #ax0.errorbar(stdp1binned[:, 0], stdp1binned[:, 1] * 100., yerr=stdp1binned[:, 2], fmt='o-',lw=2, markeredgecolor='C0')
    ax0.plot(synChange[:, 0] * 1000., synChange[:, 1] * 100.,lw=2,c='C0',label='fit: reg. pairs')
    #ax0.errorbar(stdpIrr1HzBinned[:, 0] * 1000., stdpIrr1HzBinned[:, 2] * 100., xerr=stdpIrr1HzBinned[:, 1] * 1000., yerr=stdpIrr1HzBinned[:, 3] * 100., fmt='o-',c='C3',lw=1.5,label='binned exp. data')
    if irrData:
        ax0.plot(modelNew[:, 0] * 1000., modelNew[:, 10] * 100.,lw=2,c='C2', label='prediction: irreg. pairs')
        ax0.plot(modelBurstsNew[:, 0] * 1000., modelBurstsNew[:, 10] * 100.,lw=2,c='C4', label='prediction: bursts')
        ax0.plot(modelIndNew[:, 0] * 1000., modelIndNew[:, 10] * 100.,lw=2,c='C5', label='prediction: ind. spikes')

    #for i in range(len(allData)):
    #    if allData[i][1] == 1:
    #        ax0.plot(allData[i][3] * 1000., allData[i][6] * 100., 'o', ms=4, c='0.5', label='exp. data' if i == 0 else None)

        # ax0.plot(modelNewReg[:, 0] * 1000., modelNewReg[:, 10] * 100. / 0.5,ls='--')
    # if oldNew == 'old':
    #    #print len(-model1Hz[:,3]+2.*synapticChange.D*1000.), -model1Hz[:,3]+2.*synapticChange.D*1000.
    #    ax0.plot(-model1Hz[:,3]+2.*synapticChange.D*1000.,model1Hz[:,18]*100.)
    #
    #    mask = ((-model1Hz[:,3]+2.*synapticChange.D*1000.)>15.) & ((-model1Hz[:,3]+2.*synapticChange.D*1000.)<20.)
    #    amountOfChange[dsN,2] = mean(model1Hz[mask,18]*100.)

    # else:
    # print len(modelNew[:,0]*1000.), modelNew[:,0]*1000.
    # ax0.plot(modelNew[:, 0] * 1000., modelNew[:, 1] * 100.)
    # ax0.plot(-50,modelNewFreq[4,5]*100.,'o')

    # mask = (modelNew[:, 0] * 1000. > 15.) & (modelNew[:, 0] * 1000. < 20.)
    # amountOfChange[dsN, 2] = mean(modelNew[mask, 1] * 100.)

    # removes upper and right axes
    # and moves left and bottom axes away
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_position(('outward', 10))
    ax0.spines['left'].set_position(('outward', 10))
    ax0.yaxis.set_ticks_position('left')
    ax0.xaxis.set_ticks_position('bottom')

    ax0.set_xlim(-240, 240)
    ax0.set_ylim(60, 200)
    # legends and labels
    plt.legend(loc=(0.55,0.65), frameon=False)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)

    plt.xlabel(r'time lag $\Delta t$ (ms)')
    plt.ylabel('change in synaptic strength')


    # panel 1 #############################################
    ax1 = plt.subplot(gs[1])

    # title
    ax1.set_title('2.5 - 3 Hz')

    # diplay of data
    ax1.axhline(y=100, ls='--', color='0.7', lw=2)
    ax1.axvline(x=0, ls='--', color='0.7', lw=2)
    # ax1.axvline(x=-100,ls='--',color='turquoise',lw=2)
    # ax1.axvline(x=100,ls='--',color='turquoise',lw=2)
    # ax1.axvline(x=-250,ls='--',color='turquoise',lw=2)
    # ax1.axvline(x=250,ls='--',color='turquoise',lw=2)
    #ax1.plot(stdp3Hz[:, 0], stdp3Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
    #ax1.errorbar(stdp3binned[:, 0], stdp3binned[:, 1] * 100., yerr=stdp3binned[:, 2], fmt='o-',lw=2, markeredgecolor='C0', label=r'exp. data: reg. pairs')
    ax1.plot(synChange[:, 0] * 1000., synChange[:, 2] * 100., c='C0', label=r'model fit: reg. pairs',lw=2)
    #ax1.errorbar(stdpIrr3HzBinned[:, 0] * 1000., stdpIrr3HzBinned[:, 2] * 100., xerr=stdpIrr3HzBinned[:, 1] * 1000., yerr=stdpIrr3HzBinned[:, 3] * 100., fmt='o-',
    #             c='C3',label=r'exp. data: irr. pairs')
    if irrData:
        ax1.plot(modelNew[:, 0] * 1000., modelNew[:, 11] * 100., lw=2,c='C2',label=r'model fit: irregular pairs')
        ax1.plot(modelBurstsNew[:, 0] * 1000., modelBurstsNew[:, 11] * 100.,c='C4', label=r'model fit: irregular bursts')
        ax1.plot(modelIndNew[:, 0] * 1000., modelIndNew[:, 11] * 100.,c='C5', label=r'model fit: irregular ind. spikes')
    #for i in range(len(allData)):
    #    if allData[i][1] == 3:
    #        ax1.plot(allData[i][3] * 1000., allData[i][6] * 100., 'o', ms=4, c='0.5', label='exp. data: irregular pairs' if i == 0 else None)

        # ax1.plot(modelNewReg[:, 0] * 1000., modelNewReg[:, 11] * 100. / 0.5,ls='--')  # if 'old'==oldNew:
    # ax1.plot(-model3Hz[:,3]+2.*synapticChange.D*1000.,model3Hz[:,18]*100.)
    # else:
    # ax1.plot(modelNew[:, 0] * 1000., modelNew[:, 2] * 100., label=r'model prediction: irr. pairs, Ca$_{\rm ext} = 2$ mM')

    # removes upper and right axes
    # and moves left and bottom axes away
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.spines['left'].set_position(('outward', 10))
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    #ax1.set_ylim(ymax=350)
    ax1.set_xlim(-260, 260)
    ax1.set_ylim(20, 320)
    # legends and labels
    #plt.legend(frameon=False, loc=1)
    #leg = plt.gca().get_legend()
    #ltext = leg.get_texts()
    #plt.setp(ltext, fontsize=8)

    # plt.legend(loc=1,frameon=False)

    plt.xlabel(r'time lag $\Delta t$ (ms)')
    # plt.ylabel('change in synaptic strength')


    fname = 'fig_modelPredictionBurstsIndividual_%s' % paraName  # os.path.basename(__file__)

    savefig(figDir + fname + '.png')
    savefig(figDir + fname + '.pdf')  # clf()  # dsN += 1


##############################################################################
# instance of synaptic Change and figure class  

dataDir = 'experimental_data/'
figDir = 'publicationFigures/'

generateFitRegDataFig('VenancesBin0',dataDir, figDir, modelV='bin')
