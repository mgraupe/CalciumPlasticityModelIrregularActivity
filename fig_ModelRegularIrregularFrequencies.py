from pylab import *
import os
import pickle

from synUtils import *
import params as par
import parameter_fit_solutions as pfs


################################################################################################
def generateFitRegDataFig(paraName, dataDir, figDir,modelV=None):

    w0 = 0.5
    exec('paraOpt = pfs.%s' % paraName)
    # synapticChange.choseParameterSet(analyticalLocation + 'parameters.par', fromFile=True)
    stdp1Hz = np.loadtxt(dataDir + 'STDP_1Hz_100pairings.dat')
    stdp3Hz = np.loadtxt(dataDir + 'STDP_2.5-3Hz_100pairings.dat')
    stdp5Hz = np.loadtxt(dataDir + 'STDP_5Hz_100pairings.dat')
    stdp10Hz = np.loadtxt(dataDir + 'STDP_10Hz_100pairings.dat')

    allData = pickle.load(open(dataDir + "allIrregularData.p", "rb"))

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
    gs.update(wspace=0.2, hspace=0.35)

    #fig.suptitle(r'STDP data, %s, $C_{\rm pre} = %s$, $C_{\rm post} = %s$, (RMS = %s)' % (self.modelVersion, np.round(sChange.Cpre, 4), np.round(sChange.Cpost, 4), np.round(paraOpt[1], 4)),
    #    fontsize=14)
    # possibly change outer margins of the figure
    plt.subplots_adjust(left=0.1, right=0.96, top=0.92, bottom=0.16)

    # sub-panel enumerations
    plt.figtext(0.01, 0.93, 'A',clip_on=False,color='black', weight='bold',size=22)
    plt.figtext(0.5, 0.93, 'B',clip_on=False,color='black', weight='bold',size=22)
    # plt.figtext(0.06, 0.47, 'C',clip_on=False,color='black', weight='bold',size=22)
    # plt.figtext(0.47, 0.47, 'D',clip_on=False,color='black', weight='bold',size=22)

    lineStyles = ['-',':','-.','--']
    # third sub-plot #######################################################
    ax4 = plt.subplot(gs[0])

    # title
    #ax4.set_title('regular : all frequency solutions')

    # diplay of data
    ax4.axhline(y=100, ls='--', color='0.7', lw=2)
    ax4.axvline(x=0, ls='--', color='0.7', lw=2)
    # ax4.plot(stdp5Hz[:, 0], stdp5Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
    # ax4.errorbar(stdp5binned[:, 0], stdp5binned[:, 1], yerr=stdp5binned[:, 2], fmt='o-', markeredgecolor='C0')
    for i in range(len(sChange.stimFrequencies)):
        ax4.plot(synChange[:, 0] * 1000., synChange[:, i + 1] * 100.,lw=2,alpha=(1-i*0.3),c='C0', label=str(sChange.stimFrequencies[i]))
    # ax4.plot(synChange[:, 0] * 1000., synChange[:, 2] * 100.)
    # ax4.plot(synChange[:, 0] * 1000., synChange[:, 3] * 100.)

    # if 'old'==oldNew:
    # ax2.plot(-model5Hz[:,3]+2.*synapticChange.D*1000.,model5Hz[:,18]*100.)
    # else:
    # ax2.plot(modelNew[:, 0] * 1000., modelNew[:, 3] * 100.)

    # removes upper and right axes
    # and moves left and bottom axes away
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_position(('outward', 10))
    ax4.spines['left'].set_position(('outward', 10))
    ax4.yaxis.set_ticks_position('left')
    ax4.xaxis.set_ticks_position('bottom')

    ax4.set_xlim(-260, 260)
    ax4.set_ylim(75,185)
    # legends and labels
    plt.legend(loc=(0.82,0.4), frameon=True)

    plt.xlabel(r'time lag $\Delta t$ (ms)')
    plt.ylabel('change in synaptic strength')

    # third sub-plot #######################################################
    if irrData:
        ax5 = plt.subplot(gs[1])

        # title
        #ax5.set_title('irregular : all frequency solutions')

        # diplay of data
        ax5.axhline(y=100, ls='--', color='0.7', lw=2)
        ax5.axvline(x=0, ls='--', color='0.7', lw=2)
        # ax4.plot(stdp5Hz[:, 0], stdp5Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
        # ax4.errorbar(stdp5binned[:, 0], stdp5binned[:, 1], yerr=stdp5binned[:, 2], fmt='o-', markeredgecolor='C0')
        for i in range(len(sChange.stimFrequencies)):
            ax5.plot(modelNew[:, 0] * 1000., modelNew[:, (10 + i)] * 100.,alpha=(1-i*0.3),c='C2',lw=2,label=str(sChange.stimFrequencies[i]))  # ax5.plot(synChange[:, 0] * 1000., synChange[:, i + 1] * 100., label=str(self.stimFrequencies[i]))
        # ax4.plot(synChange[:, 0] * 1000., synChange[:, 2] * 100.)
        # ax4.plot(synChange[:, 0] * 1000., synChange[:, 3] * 100.)

        # if 'old'==oldNew:
        # ax2.plot(-model5Hz[:,3]+2.*synapticChange.D*1000.,model5Hz[:,18]*100.)
        # else:
        # ax2.plot(modelNew[:, 0] * 1000., modelNew[:, 3] * 100.)

        # removes upper and right axes
        # and moves left and bottom axes away
        ax5.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['bottom'].set_position(('outward', 10))
        ax5.spines['left'].set_position(('outward', 10))
        ax5.yaxis.set_ticks_position('left')
        ax5.xaxis.set_ticks_position('bottom')

        ax5.set_xlim(-260, 260)
        ax5.set_ylim(75,185)
        # legends and labels
        plt.legend(loc=(0.82,0.4), frameon=True)

        plt.xlabel(r'time lag $\Delta t$ (ms)')  # plt.ylabel('change in synaptic strength')

    ## save figure ############################################################
    # fname = 'regularDataFitVenance' #_dataSet#' + str(dataSetNumber)


    fname = 'fig_regularIrregularFrequencies_%s' % paraName

    savefig(figDir + fname + '.png')
    savefig(figDir + fname + '.pdf')  # clf()  # dsN += 1



##############################################################################
# instance of synaptic Change and figure class


dataDir = 'experimental_data/'
figDir = 'publicationFigures/'


generateFitRegDataFig('VenancesBin0',dataDir, figDir, modelV='bin')



