from pylab import *
import os
import pickle
import pdb

from synUtils import *
import params as par
import parameter_fit_solutions as pfs


################################################################################################
def generateFitRegDataFig(paraName, allData, figDir,burstData,distributionBursts, modelV=None):

    w0 = 0.5
    exec('paraOpt = pfs.%s' % paraName)
    # synapticChange.choseParameterSet(analyticalLocation + 'parameters.par', fromFile=True)
    #stdp1Hz = np.loadtxt(dataDir + 'STDP_1Hz_100pairings.dat')
    #stdp3Hz = np.loadtxt(dataDir + 'STDP_2.5-3Hz_100pairings.dat')
    #stdp5Hz = np.loadtxt(dataDir + 'STDP_5Hz_100pairings.dat')
    #stdp10Hz = np.loadtxt(dataDir + 'STDP_10Hz_100pairings.dat')

    #allData = pickle.load(open(dataDir + "allIrregularData.p", "rb")) #, encoding='latin1')

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

    #######################################################
    # plot data
    fig_width = 10  # width in inches
    fig_height = 7  # height in inches
    fig_size = [fig_width, fig_height]
    params = {'axes.labelsize': 11, 'axes.titlesize': 11, 'font.size': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11, 'figure.figsize': fig_size,  # 'savefig.dpi': 600,
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
    gs = gridspec.GridSpec(2, 1  # width_ratios=[1,1.2],
                           # height_ratios=[1,1]
                           )

    # define vertical and horizontal spacing between panels
    gs.update(wspace=0.3, hspace=0.35)

    #fig.suptitle(r'STDP data, %s, $C_{\rm pre} = %s$, $C_{\rm post} = %s$, (RMS = %s)' % (modelVersion, np.round(sChange.Cpre, 4), np.round(sChange.Cpost, 4), np.round(paraOpt[1], 4)),
    #    fontsize=14)
    # possibly change outer margins of the figure
    plt.subplots_adjust(left=0.1, right=0.96, top=0.92, bottom=0.1)

    # sub-panel enumerations
    plt.figtext(0.01, 0.95, 'L',clip_on=False,color='black', weight='bold',size=12)
    plt.figtext(0.53, 0.95, 'M',clip_on=False,color='black', weight='bold',size=12)
    plt.figtext(0.01, 0.455, 'C',clip_on=False,color='black', weight='bold',size=12)
    plt.figtext(0.32, 0.455, 'D',clip_on=False,color='black', weight='bold',size=12)
    plt.figtext(0.6, 0.455, 'E',clip_on=False,color='black', weight='bold',size=12)

    # panel 0 #######################################################
    gssub0 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], wspace=0.3)
    ax0 = plt.subplot(gssub0[0])

    # title
    ax0.set_title('1 spk/s')

    # diplay of data
    ax0.axhline(y=100, ls='--', color='0.7', lw=2)
    ax0.axvline(x=0, ls='--', color='0.7', lw=2)
    # ax0.axvline(x=-50,ls='--',color='turquoise',lw=2)
    # ax0.axvline(x=25,ls='--',color='turquoise',lw=2)
    # ax0.axvline(x=-200,ls='--',color='turquoise',lw=2)
    # ax0.axvline(x=200,ls='--',color='turquoise',lw=2)
    #ax0.plot(stdp1Hz[:, 0], stdp1Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
    #ax0.errorbar(stdp1binned[:, 0], stdp1binned[:, 1] * 100., yerr=stdp1binned[:, 2], fmt='o-',lw=2, markeredgecolor='C0')
    ax0.plot(allData[6][2][:, 0]*1000., allData[6][2][:, 1] * 100., 'o', ms=7,alpha=0.5, c='white',markeredgewidth=1,markeredgecolor='#7575ff', label='burst exp. data',zorder=10)
    ax0.errorbar(allData[4][3][:, 0] * 1000., allData[4][3][:, 1] * 100., xerr=allData[4][3][:, 2] * 1000., yerr=allData[4][3][:, 3] * 100., alpha=0.3, c='C1', lw=2, fmt='o-',label='binned irregular exp. data')
    ax0.errorbar(allData[6][3][:, 0]*1000., allData[6][3][:, 1] * 100., xerr=allData[6][3][:, 2]*1000., yerr=allData[6][3][:, 3] * 100., c='#7575ff', lw=2, fmt='o-', zorder=20, label='binned burst exp. data')

    print('1 spk/s burst data (binned):')
    print(allData[6][3])
    #ax0.plot(synChange[:, 0] * 1000., synChange[:, 1] * 100.,lw=2.5,c='C2',label='fit: reg. pairs')
    #ax0.plot(modelNew[:, 0] * 1000., modelNew[:, 10] * 100.,lw=2.5,c='C5', label='prediction: irreg. pairs')
    #ax0.plot(modelBurstsNew[:, 0] * 1000., modelBurstsNew[:, 10] * 100.,lw=2.5,c='C4', label='prediction: bursts')

    #ax0.plot(modelBurstsNew[:, 0] * 1000., modelBurstsNew[:, 10] * 100., label='bursts')
    #ax0.plot(modelIndNew[:, 0] * 1000., modelIndNew[:, 10] * 100., label='ind. spikes')

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

    ax0.set_xlim(-220, 240)
    ax0.set_ylim(60, 240)
    # legends and labels
    plt.legend(loc=(0.55,0.5), frameon=False)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)

    plt.xlabel(r'time lag $\Delta t$ (ms)')
    plt.ylabel('change in synaptic strength')

    # panel 1 #############################################
    ax1 = plt.subplot(gssub0[1])

    # title
    ax1.set_title('3 spk/s')

    # diplay of data
    ax1.axhline(y=100, ls='--', color='0.7', lw=2)
    ax1.axvline(x=0, ls='--', color='0.7', lw=2)
    # ax1.axvline(x=-100,ls='--',color='turquoise',lw=2)
    # ax1.axvline(x=100,ls='--',color='turquoise',lw=2)
    # ax1.axvline(x=-250,ls='--',color='turquoise',lw=2)
    # ax1.axvline(x=250,ls='--',color='turquoise',lw=2)
    #ax1.plot(stdp3Hz[:, 0], stdp3Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
    #ax1.errorbar(stdp3binned[:, 0], stdp3binned[:, 1] * 100., yerr=stdp3binned[:, 2], fmt='o-',lw=2, markeredgecolor='C0', label=r'exp. data: reg. pairs')
    ax1.plot(allData[7][2][:, 0]*1000., allData[7][2][:, 1] * 100., 'o', ms=5,alpha=0.5, c='white',markeredgewidth=1,markeredgecolor='#7575ff', label='burst exp. data',zorder=10)
    ax1.errorbar(allData[5][3][:, 0] * 1000., allData[5][3][:, 1] * 100., xerr=allData[5][3][:, 2] * 1000., yerr=allData[5][3][:, 3] * 100., c='C1', lw=2, alpha=0.3,fmt='o-',label='binned irregular exp. data')
    ax1.errorbar(allData[7][3][:, 0]*1000., allData[7][3][:, 1] * 100., xerr=allData[7][3][:, 2]*1000., yerr=allData[7][3][:, 3] * 100., c='#7575ff', lw=2, fmt='o-',zorder=20, label='binned burst exp. data')

    print('3 spk/s burst data (binned):')
    print(allData[7][3])

    #ax1.plot(synChange[:, 0] * 1000., synChange[:, 2] * 100., c='C2',lw=2.5, label=r'model fit: reg. pairs')
    #ax1.plot(modelNew[:, 0] * 1000., modelNew[:, 11] * 100., lw=2.5,c='C5',label=r'model fit: irregular pairs')
    #ax1.plot(modelBurstsNew[:, 0] * 1000., modelBurstsNew[:, 11] * 100.,lw=2.5,c='C4', label=r'model fit: irregular bursts')


    # removes upper and right axes
    # and moves left and bottom axes away
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('outward', 10))
    ax1.spines['left'].set_position(('outward', 10))
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    ax1.set_xlim(-220, 240)
    ax1.set_ylim(60, 240)
    # legends and labels
    #plt.legend(frameon=False, loc=1)
    #leg = plt.gca().get_legend()
    #ltext = leg.get_texts()
    #plt.setp(ltext, fontsize=8)

    # plt.legend(loc=1,frameon=False)

    plt.xlabel(r'time lag $\Delta t$ (ms)')
    # plt.ylabel('change in synaptic strength')

    # panel 0 #######################################################
    gssub1 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.45, width_ratios=[1,1,1.5])
    ax2 = plt.subplot(gssub1[0])

    # title
    #ax2.set_title('1 Hz')
    pairRate = [2.5,2.,3.,4.]
    idx0 = 0
    idx1 = 2
    #pdb.set_trace()
    for n in range(len(allData[6][4])):
        #print('number of pre bursts : ', allData[6][4][n][7][0][2]) #[6][2])
        nPreBursts =  allData[6][4][n][7][0][2]
        #print('pre-spikes in bursts : ',allData[6][4][n][7][0][3]) #[6][3])
        preSpikesInBursts = np.sum(allData[6][4][n][7][0][3])
        #print('number of post bursts : ', allData[6][4][n][7][1][2])  # [6][2])
        nPostBursts = allData[6][4][n][7][1][2]
        #print('post-spikes in bursts : ', allData[6][4][n][7][1][3])  # [6][3])
        postSpikesInBursts = np.sum(allData[6][4][n][7][1][3])
        ax2.plot([2.4,2.6],[preSpikesInBursts,postSpikesInBursts],'o-',ms=2,alpha=0.6,c='0.5')

    for n in range(len(allData[7][4])):
        #print('number of pre bursts : ', allData[6][4][n][7][0][2]) #[6][2])
        nPreBursts =  allData[7][4][n][7][0][2]
        #print('pre-spikes in bursts : ',allData[6][4][n][7][0][3]) #[6][3])
        preSpikesInBursts = np.sum(allData[7][4][n][7][0][3])
        #print('number of post bursts : ', allData[6][4][n][7][1][2])  # [6][2])
        nPostBursts = allData[7][4][n][7][1][2]
        #print('post-spikes in bursts : ', allData[6][4][n][7][1][3])  # [6][3])
        postSpikesInBursts = np.sum(allData[7][4][n][7][1][3])
        ax2.plot([2.9,3.1],[preSpikesInBursts,postSpikesInBursts],'o-',ms=2,alpha=0.6,c='0.5')
        if preSpikesInBursts < 20.:
            print(allData[7][4][n])

    ax2.errorbar(pairRate[idx0], np.mean(burstData[idx0, 0, :]),yerr=np.std(burstData[idx0, 0, :]),marker='o', color='C0')
    ax2.errorbar(pairRate[idx1],np.mean(burstData[idx1, 0, :]),yerr=np.std(burstData[idx1, 0, :]),marker='o', color='C0')


    # removes upper and right axes
    # and moves left and bottom axes away
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_position(('outward', 10))
    ax2.spines['left'].set_position(('outward', 10))
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')

    #ax2.set_xlim(-260, 260)
    #ax2.set_ylim(55, 320)
    # legends and labels
    #majorLocator_x = MultipleLocator(2)
    #ax2.xaxis.set_major_locator(majorLocator_x)
    ax2.xaxis.set_ticks([2.5,3])
    ax2.xaxis.set_ticklabels(['1','3'])
    #plt.legend(loc=(0.55, 0.5), frameon=False)
    #leg = plt.gca().get_legend()
    #ltext = leg.get_texts()
    #plt.setp(ltext, fontsize=8)

    plt.xlabel(r'rate (spk/s)')
    plt.ylabel('total number of spikes in bursts')
    ########################################################################
    ax2 = plt.subplot(gssub1[1])

    # title
    # ax2.set_title('1 Hz')
    pairRate = [2.5, 2., 3., 4.]
    idx0 = 0
    idx1 = 2
    # pdb.set_trace()
    for n in range(len(allData[6][4])):
        # print('number of pre bursts : ', allData[6][4][n][7][0][2]) #[6][2])
        nPreBursts = allData[6][4][n][7][0][2]
        # print('pre-spikes in bursts : ',allData[6][4][n][7][0][3]) #[6][3])
        preSpikesInBursts = np.sum(allData[6][4][n][7][0][3])
        # print('number of post bursts : ', allData[6][4][n][7][1][2])  # [6][2])
        nPostBursts = allData[6][4][n][7][1][2]
        # print('post-spikes in bursts : ', allData[6][4][n][7][1][3])  # [6][3])
        postSpikesInBursts = np.sum(allData[6][4][n][7][1][3])
        ax2.plot([2.4, 2.6], [nPreBursts, nPostBursts], 'o-', ms=2, alpha=0.6, c='0.5')

    for n in range(len(allData[7][4])):
        # print('number of pre bursts : ', allData[6][4][n][7][0][2]) #[6][2])
        nPreBursts = allData[7][4][n][7][0][2]
        # print('pre-spikes in bursts : ',allData[6][4][n][7][0][3]) #[6][3])
        preSpikesInBursts = np.sum(allData[7][4][n][7][0][3])
        # print('number of post bursts : ', allData[6][4][n][7][1][2])  # [6][2])
        nPostBursts = allData[7][4][n][7][1][2]
        # print('post-spikes in bursts : ', allData[6][4][n][7][1][3])  # [6][3])
        postSpikesInBursts = np.sum(allData[7][4][n][7][1][3])
        ax2.plot([2.9, 3.1], [nPreBursts, nPostBursts], 'o-', ms=2, alpha=0.6, c='0.5')
        #if preSpikesInBursts < 20.:
        #    print(allData[7][4][n])

    ax2.errorbar(pairRate[idx0], np.mean(burstData[idx0, 2, :]), yerr=np.std(burstData[idx0, 2, :]), marker='o', color='C0')
    ax2.errorbar(pairRate[idx1], np.mean(burstData[idx1, 2, :]), yerr=np.std(burstData[idx1, 2, :]), marker='o', color='C0')

    # removes upper and right axes
    # and moves left and bottom axes away
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_position(('outward', 10))
    ax2.spines['left'].set_position(('outward', 10))
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')

    # ax2.set_xlim(-260, 260)
    # ax2.set_ylim(55, 320)
    # legends and labels
    # majorLocator_x = MultipleLocator(2)
    # ax2.xaxis.set_major_locator(majorLocator_x)
    ax2.xaxis.set_ticks([2.5, 3])
    ax2.xaxis.set_ticklabels(['1','3'])
    # plt.legend(loc=(0.55, 0.5), frameon=False)
    # leg = plt.gca().get_legend()
    # ltext = leg.get_texts()
    # plt.setp(ltext, fontsize=8)

    plt.xlabel(r'rate (spk/s)')
    plt.ylabel('total number bursts')
    ########################################################################
    ax5 = plt.subplot(gssub1[2])
    # title
    # ax1.set_title('3 Hz')
    labels = ['1 spk/s','3 spk/s']
    x = np.arange(len(labels))
    width = 0.12
    spacing = 0.3
    barLocations1 = x[0] + np.arange(-2,3,1)*spacing/2
    barLocations3 = x[1] + np.arange(-2,3,1)*spacing/2
    for n in range(len(allData[6][4])):
        # print('number of pre bursts : ', allData[6][4][n][7][0][2]) #[6][2])
        #nPreBursts = allData[6][4][n][7][0][2]
        # print('pre-spikes in bursts : ',allData[6][4][n][7][0][3]) #[6][3])
        preBins = np.bincount(allData[6][4][n][7][0][3])
        lPre = len(preBins)-2
        #pdb.set_trace()
        # print('number of post bursts : ', allData[6][4][n][7][1][2])  # [6][2])
        #nPostBursts = allData[6][4][n][7][1][2]
        # print('post-spikes in bursts : ', allData[6][4][n][7][1][3])  # [6][3])
        postBins = np.bincount(allData[6][4][n][7][1][3])
        lPost = len(postBins)-2
        ax5.plot(barLocations1[:lPre], preBins[2:], 'o-', ms=2, alpha=0.4, c='black')
        ax5.plot(barLocations1[:lPost], postBins[2:], 'o-', ms=2, alpha=0.4, c='turquoise')

    for n in range(len(allData[7][4])):
        # print('number of pre bursts : ', allData[6][4][n][7][0][2]) #[6][2])
        #nPreBursts = allData[6][4][n][7][0][2]
        # print('pre-spikes in bursts : ',allData[6][4][n][7][0][3]) #[6][3])
        preBins = np.bincount(allData[7][4][n][7][0][3])
        lPre = len(preBins[:7])-2
        #pdb.set_trace()
        # print('number of post bursts : ', allData[6][4][n][7][1][2])  # [6][2])
        #nPostBursts = allData[6][4][n][7][1][2]
        # print('post-spikes in bursts : ', allData[6][4][n][7][1][3])  # [6][3])
        postBins = np.bincount(allData[7][4][n][7][1][3])
        lPost = len(postBins[:7])-2
        ax5.plot(barLocations3[:lPre], preBins[2:7], 'o-', ms=2, alpha=0.4, c='black')
        ax5.plot(barLocations3[:lPost], postBins[2:7], 'o-', ms=2, alpha=0.4, c='turquoise')


    rects1 = ax5.bar(x - 2 * spacing / 2, distributionBursts[[0,2], 2], width, label='2')
    rects2 = ax5.bar(x - spacing / 2, distributionBursts[[0,2], 3], width, label='3')
    rects3 = ax5.bar(x, distributionBursts[[0,2], 4], width, label='4')
    rects4 = ax5.bar(x + spacing / 2, distributionBursts[[0,2], 5], width, label='5')
    rects5 = ax5.bar(x + 2 * spacing / 2, distributionBursts[[0,2], 6], width, label='6')
    # rects5 = ax5.bar(x + 2*spacing/2, distributionBursts[:,6], width, label='6')

    # removes upper and right axes
    # and moves left and bottom axes away
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_position(('outward', 10))
    ax5.spines['left'].set_position(('outward', 10))
    ax5.yaxis.set_ticks_position('left')
    ax5.xaxis.set_ticks_position('bottom')

    ax5.set_xticks(x)
    ax5.set_xticklabels(labels)

    ax5.set_ylabel('occurrence of bursts with X spikes')

    plt.legend(frameon=False, loc=(0.85, 0.6))
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)

    ########################################################################
    fname = 'fig_modelAndDataBurstStim_%s' % paraName  # os.path.basename(__file__)

    savefig(figDir + fname + '.png')
    savefig(figDir + fname + '.pdf')  # clf()  # dsN += 1


##############################################################################
# instance of synaptic Change and figure class  
dataDir = '/media/HDnyc_data/data_analysis/SinglePlasticityTraces/'
figDir = 'publicationFigures/'

experimentalPlasticityData = pickle.load(open(dataDir+'experimentalPlasticityData.p', 'rb'))

burstData = np.load(open('burstData.npy','rb'))
distributionBursts = np.load(open('distributionBursts.npy','rb'))
#pdb.set_trace()
generateFitRegDataFig('VenancesBin0',experimentalPlasticityData, figDir,burstData,distributionBursts, modelV='bin')
