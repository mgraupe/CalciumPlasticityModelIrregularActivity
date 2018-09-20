from pylab import *
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import time

from timeAboveThreshold import timeAboveThreshold
from synapticChange import synapticChange
import parameter_fit_solutions as pfs

class synUtils(): # synUtils(thetaD,thetaP,nonlinear,Npresentations,w0)
    ###############################################################################
    def __init__(self, dataSet='venance',modelV='multiplicative'):
        #self.thetaD  = thetaD
        #self.thetaP  = thetaP
        # determine eta based on nonlinearity factor and amplitudes
        #self.nonlinear = nonlinear
        #self.Npresentations = Npresentations
        #self.w0 = w0
        self.modelVersion = modelV
        #self.nonlinear = 1.
        # read in experimental data
        #dataDir = '/home/mgraupe/theobio/network_1/fit_all_models/calcium_nonlinear_python_parameter_search/experimental_data/'
        self.dataDir = 'experimental_data/'
        self.figDir = 'FigsSimResults/'

        if dataSet == 'sjoestroem':
            jesperReg = loadtxt(self.dataDir+'sjoestroem_regular_all.dat')
            jesperStoch = loadtxt(self.dataDir+'sjoestroem_stochastic.dat')

            self.xDataReg = jesperReg[:,[0,1]]
            self.xDataReg[:,1] = self.xDataReg[:,1]/1000. # everything in sec
            self.yDataReg = jesperReg[:,2]+1. # Sjoestroem's data is normalized to 0
            self.sigmaDataReg = jesperReg[:,3]


            self.xDataStoch = jesperStoch[:,0]
            self.yDataStoch = jesperStoch[:,1]+1. # Sjoestroem's data is normalized to 0
            self.sigmaDataStoch = jesperStoch[:,2]
        elif dataSet == 'venance':
            self.stimFrequencies = [1,3,5,10]
            #self.fitWeights = [4.,1.,1.,1.]
            self.Npresentations = 100
            self.rawData1Hz  = loadtxt(self.dataDir+'STDP_1Hz_100pairings_binned.dat')
            self.rawData3Hz = loadtxt(self.dataDir+'STDP_2.5-3Hz_100pairings_binned.dat')
            self.rawData5Hz = loadtxt(self.dataDir+'STDP_5Hz_100pairings_binned.dat')
            self.rawData10Hz = loadtxt(self.dataDir+'STDP_10Hz_100pairings_binned.dat')

            #self.fitData = {}
            #self.fitData[0] = {'freq':1.,'data':rawData1Hz}
            #self.fitData[1] = {'freq':3.,'data':rawData3Hz}
            #self.fitData[2] = {'freq':5.,'data':rawData5Hz}
            #self.fitData[3] = {'freq':10.,'data':rawData10Hz}

    ################################################################################################
    def generateFig(self, paraOpt, figName = None):
        
        ####################################
        # calculate solution
        freq = linspace(0.1,50.,500)
        synChange = zeros((len(freq),2))
        synChangeStoch = zeros((len(freq)))

        deltaTs = linspace(-0.05,0.05,1001)
        synChange2 = zeros(len(deltaTs))

        for i in range(len(freq)):
            #calculateChangeInSynapticStrength(frequency,deltaT,params):
            synChange[i,0] = self.calculateChangeInSynapticStrength(freq[i],0.01,paraOpt[0])
            synChange[i,1] = self.calculateChangeInSynapticStrength(freq[i],-0.01,paraOpt[0])
            #synChangeStoch[i] = self.calculateChangeInSynapticStrengthStochastic(freq[i],paraOpt[0],[-0.015,0.015])
            synChangeStoch[i] = self.calculateChangeInSynapticStrengthStochastic(freq[i],paraOpt[0],[-0.015,0.015])
            
        for i in range(len(deltaTs)):
            synChange2[i] = self.calculateChangeInSynapticStrength(0.1,deltaTs[i],paraOpt[0])


        fig_width = 5 # width in inches
        fig_height = 12  # height in inches
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
        gs = gridspec.GridSpec(3, 1,
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
        ax0.set_title('regular Sjoestroem \n tauCa=%.6f, Cpre=%.6f, Cpost=%.6f\n gammaD=%.6f, gammaP=%.6f, tau=%.6f, D=%.6f\n chi2=%.3f' % (paraOpt[0][0],paraOpt[0][1],paraOpt[0][2],paraOpt[0][3],paraOpt[0][4],paraOpt[0][5],paraOpt[0][6],paraOpt[1]),y=1.05,fontsize=11)
        #, chi2='+str(paraOpt[1]))

        # diplay of data
        ax0.axhline(y=1.,c='0.7')
        ax0.errorbar(self.xDataReg[:,0][::2],self.yDataReg[::2],yerr=self.sigmaDataReg[::2],fmt='s',color='red',clip_on=False)
        ax0.errorbar(self.xDataReg[:,0][1::2],self.yDataReg[1::2],yerr=self.sigmaDataReg[1::2],fmt='o',color='blue',clip_on=False)
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
        ax0 = plt.subplot(gs[1])

        # title
        ax0.set_title('stochastic Sjoestroem')

        # diplay of data
        ax0.axhline(y=1.,c='0.7')
        ax0.errorbar(self.xDataStoch,self.yDataStoch,yerr=self.sigmaDataStoch,fmt='s',color='green',clip_on=False)
        ax0.plot(freq,synChangeStoch,color='green')

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
        ax1 = plt.subplot(gs[2])

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
        ## save figure ############################################################
        if figName:
            fname = 'fitFigures/' + figName #os.path.basename(__file__)
        else:
            fname = 'synChangeCa-ModelFit' #os.path.basename(__file__)

        savefig(self.figDir + fname+'.png')
        savefig(self.figDir + fname+'.pdf')
        # close figures to avoid out of memory
        if figName:
            plt.close(fig)
            #fig.clf()

    ################################################################################################
    def generateVenFig(self, paraName, figName=None, modelV = None):

        self.w0 = 0.5
        exec('paraOpt = pfs.%s' % paraName)
        #synapticChange.choseParameterSet(analyticalLocation + 'parameters.par', fromFile=True)
        stdp1Hz = np.loadtxt(self.dataDir + 'STDP_1Hz_100pairings.dat')
        stdp3Hz = np.loadtxt(self.dataDir + 'STDP_2.5-3Hz_100pairings.dat')
        stdp5Hz = np.loadtxt(self.dataDir + 'STDP_5Hz_100pairings.dat')
        stdp10Hz = np.loadtxt(self.dataDir+ 'STDP_10Hz_100pairings.dat')

        stdp1binned = np.copy(self.rawData1Hz)
        stdp3binned = np.copy(self.rawData3Hz)
        stdp5binned = np.copy(self.rawData5Hz)
        stdp10binned = np.copy(self.rawData10Hz)
        ##########################################################
        # generate analytical results for regular pairs
        if os.path.isfile(self.figDir+'irregularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName):
            irrData = True
            modelNew = np.load(self.figDir+'irregularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName)
            modelNewReg = np.load(self.figDir+'regularSpikePairs_vs_deltaT_differentFreqs_%s.npy' % paraName)
        else:
            irrData = False
        #Npairs = 100.
        #tat = timeAboveThreshold(synapticChange.thetaD, synapticChange.thetaP, synapticChange.tauCa, synapticChange.Cpre * CaTest / Ca0, synapticChange.Cpost * CaTest / Ca0)

        ###########################################################
        # change as function of delta t for four different frequencies
        deltaTstart = -0.3
        deltaTend = 0.3
        steps = 3001.

        #stimFreq = [1,3,5,10]

        if modelV is not None:
            if modelV == 'a':
                self.modelVersion = 'additive'
            elif modelV == 'm':
                self.modelVersion = 'multiplicative'
            else:
                print 'problem in model choice'
                sys.exit(1)

        sChange = synapticChange('Venance') #,threshold=par.thetaP)
        sChange.choseParameterSet(paraName, source='fromFile')
        # initiate class to calculate fraction of time above threshold
        tat = timeAboveThreshold(sChange.tauCa, sChange.Cpre, sChange.Cpost, sChange.thetaD, sChange.thetaP)
        print 'Parameters :', sChange.tauCa, sChange.Cpre, sChange.Cpost, sChange.thetaD, sChange.thetaP

        deltaT = linspace(deltaTstart, deltaTend, steps)
        synChange = zeros((len(deltaT), len(sChange.stimFrequencies) + 1))
        for n in range(len(sChange.stimFrequencies)):
            #frequency = stimFreq[n]
            interval = 1. / sChange.stimFrequencies[n]
            for i in range(len(deltaT)):
                (alphaD, alphaP) = tat.spikePairFrequencyNonlinear(deltaT[i] - sChange.D, sChange.stimFrequencies[n])
                # print dT, preRate, alphaD, alphaP
                sChange.changeInSynapticStrength(sChange.Npresentations / sChange.stimFrequencies[n], self.w0, alphaD, alphaP)

                if n == 0:
                    synChange[i, 0] = deltaT[i]
                #  calculateChangeInSynapticStrength(self, frequency,deltaT,params):
                #sol = [0.0667179, 1.45248, 0.405039, 2.0, 15.973, 16.3457, -0.00156591]
                if self.modelVersion == 'additive':
                    synChange[i, n+1] = sChange.meanAdd/self.w0
                elif self.modelVersion == 'multiplicative':
                    synChange[i, n+1] = sChange.mean/self.w0
                #print self.stimFrequencies[n], deltaT[i], synChange[i, n+1]
            #
        #pdb.set_trace()
        #######################################################
        # plot data
        fig_width = 10  # width in inches
        fig_height = 15  # height in inches
        fig_size = [fig_width, fig_height]
        params = {'axes.labelsize': 14,
                  'axes.titlesize': 13,
                  'font.size': 11,
                  'xtick.labelsize': 11,
                  'ytick.labelsize': 11,
                  'figure.figsize': fig_size,
                  #'savefig.dpi': 600,
                  'axes.linewidth': 1.3,
                  'ytick.major.size': 4,  # major tick size in points
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
        gs = gridspec.GridSpec(3, 2  # width_ratios=[1,1.2],
                               # height_ratios=[1,1]
                               )

        # define vertical and horizontal spacing between panels
        gs.update(wspace=0.2, hspace=0.35)

        fig.suptitle(
            r'STDP data, %s, $C_{\rm pre} = %s$, $C_{\rm post} = %s$, (RMS = %s)' % (self.modelVersion,np.round(sChange.Cpre,4), np.round(sChange.Cpost,4), np.round(paraOpt[1],4)),
            fontsize=14)
        # possibly change outer margins of the figure
        plt.subplots_adjust(left=0.14, right=0.92, top=0.92, bottom=0.08)

        # sub-panel enumerations
        # plt.figtext(0.06, 0.92, 'A',clip_on=False,color='black', weight='bold',size=22)
        # plt.figtext(0.47, 0.92, 'B',clip_on=False,color='black', weight='bold',size=22)
        # plt.figtext(0.06, 0.47, 'C',clip_on=False,color='black', weight='bold',size=22)
        # plt.figtext(0.47, 0.47, 'D',clip_on=False,color='black', weight='bold',size=22)

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
        ax0.plot(stdp1Hz[:, 0], stdp1Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
        ax0.errorbar(stdp1binned[:, 0], stdp1binned[:, 1], yerr=stdp1binned[:, 2], fmt='o-', markeredgecolor='C0')
        ax0.plot(synChange[:, 0] * 1000., synChange[:, 1] * 100.)
        if irrData:
            ax0.plot(modelNew[:,0]*1000.,modelNew[:,10]*100./0.5)
            ax0.plot(modelNewReg[:, 0] * 1000., modelNewReg[:, 10] * 100. / 0.5,ls='--')
        # if oldNew == 'old':
        #    #print len(-model1Hz[:,3]+2.*synapticChange.D*1000.), -model1Hz[:,3]+2.*synapticChange.D*1000.
        #    ax0.plot(-model1Hz[:,3]+2.*synapticChange.D*1000.,model1Hz[:,18]*100.)
        #
        #    mask = ((-model1Hz[:,3]+2.*synapticChange.D*1000.)>15.) & ((-model1Hz[:,3]+2.*synapticChange.D*1000.)<20.)
        #    amountOfChange[dsN,2] = mean(model1Hz[mask,18]*100.)

        # else:
        # print len(modelNew[:,0]*1000.), modelNew[:,0]*1000.
        #ax0.plot(modelNew[:, 0] * 1000., modelNew[:, 1] * 100.)
        # ax0.plot(-50,modelNewFreq[4,5]*100.,'o')

        #mask = (modelNew[:, 0] * 1000. > 15.) & (modelNew[:, 0] * 1000. < 20.)
        #amountOfChange[dsN, 2] = mean(modelNew[mask, 1] * 100.)

        # removes upper and right axes
        # and moves left and bottom axes away
        ax0.spines['top'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.spines['bottom'].set_position(('outward', 10))
        ax0.spines['left'].set_position(('outward', 10))
        ax0.yaxis.set_ticks_position('left')
        ax0.xaxis.set_ticks_position('bottom')

        ax0.set_xlim(-300, 300)
        # legends and labels
        # plt.legend(loc=1,frameon=False)

        # plt.xlabel(r'$\Delta t$ (ms)')
        plt.ylabel('change in synaptic strength')

        ax_inset = plt.axes((0.37, 0.81, 0.13, 0.15))
        # plt.hist(residuals, fc='0.8',ec='w', lw=2)
        ax_inset.axhline(y=100, ls='--', color='0.7', lw=2)
        ax_inset.axvline(x=0, ls='--', color='0.7', lw=2)
        # ax_inset.axvline(x=-50,ls='--',color='turquoise',lw=2)
        # ax_inset.axvline(x=25,ls='--',color='turquoise',lw=2)
        ax_inset.plot(stdp1Hz[:, 0], stdp1Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
        ax_inset.errorbar(stdp1binned[:, 0], stdp1binned[:, 1], yerr=stdp1binned[:, 2], fmt='o-', markeredgecolor='C0')
        ax_inset.plot(synChange[:, 0] * 1000., synChange[:, 1] * 100.)
        if irrData:
            ax_inset.plot(modelNew[:,0]*1000.,modelNew[:,10]*100./0.5)

            # if oldNew == 'old':
        #    ax_inset.plot(-model1Hz[:,3]+2.*synapticChange.D*1000.,model1Hz[:,18]*100.)
        # else:
        #ax_inset.plot(modelNew[:, 0] * 1000., modelNew[:, 1] * 100.)
        ax_inset.set_xlim(-40, 40)

        ax_inset.spines['top'].set_visible(False)
        ax_inset.spines['right'].set_visible(False)
        ax_inset.spines['bottom'].set_position(('outward', 10))
        ax_inset.spines['left'].set_position(('outward', 10))
        ax_inset.yaxis.set_ticks_position('left')
        ax_inset.xaxis.set_ticks_position('bottom')

        ax_inset.yaxis.set_major_locator(MaxNLocator(3))
        ax_inset.xaxis.set_major_locator(MaxNLocator(3))
        # plt.xticks([-0.5, 0, 0.5],[-0.5, 0, 0.5], size=6)
        # plt.xlim((-0.5, 0.5))
        # plt.yticks([5, 10], size=6)
        # plt.xlabel("residuals", size=8)
        # ax_inset.xaxis.set_ticks_position("none")
        # ax_inset.yaxis.set_ticks_position("left")

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
        ax1.plot(stdp3Hz[:, 0], stdp3Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
        ax1.errorbar(stdp3binned[:, 0], stdp3binned[:, 1], yerr=stdp3binned[:, 2], fmt='o-', markeredgecolor='C0', label=r'exp. data: reg. pairs')
        ax1.plot(synChange[:, 0] * 1000., synChange[:, 2] * 100., label=r'model fit: reg. pairs')
        if irrData:
            ax1.plot(modelNew[:,0]*1000.,modelNew[:,11]*100./0.5, label=r'model pred.: irregular pairs')
            ax1.plot(modelNewReg[:, 0] * 1000., modelNewReg[:, 11] * 100. / 0.5,ls='--')
            # if 'old'==oldNew:
        # ax1.plot(-model3Hz[:,3]+2.*synapticChange.D*1000.,model3Hz[:,18]*100.)
        # else:
        #ax1.plot(modelNew[:, 0] * 1000., modelNew[:, 2] * 100., label=r'model prediction: irr. pairs, Ca$_{\rm ext} = 2$ mM')

        # removes upper and right axes
        # and moves left and bottom axes away
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_position(('outward', 10))
        ax1.spines['left'].set_position(('outward', 10))
        ax1.yaxis.set_ticks_position('left')
        ax1.xaxis.set_ticks_position('bottom')

        ax1.set_xlim(-300, 300)
        # legends and labels
        plt.legend(frameon=False, loc=1)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=8)

        # plt.legend(loc=1,frameon=False)

        # plt.xlabel(r'$\Delta t$ (ms)')
        # plt.ylabel('change in synaptic strength')

        # third sub-plot #######################################################
        ax2 = plt.subplot(gs[2])

        # title
        ax2.set_title('5 Hz')

        # diplay of data
        ax2.axhline(y=100, ls='--', color='0.7', lw=2)
        ax2.axvline(x=0, ls='--', color='0.7', lw=2)
        ax2.plot(stdp5Hz[:, 0], stdp5Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
        ax2.errorbar(stdp5binned[:, 0], stdp5binned[:, 1], yerr=stdp5binned[:, 2], fmt='o-', markeredgecolor='C0')
        ax2.plot(synChange[:, 0] * 1000., synChange[:, 3] * 100.)
        if irrData:
            ax2.plot(modelNew[:,0]*1000.,modelNew[:,12]*100./0.5)
            ax2.plot(modelNewReg[:, 0] * 1000., modelNewReg[:, 12] * 100. / 0.5,ls='--')
            # if 'old'==oldNew:
        # ax2.plot(-model5Hz[:,3]+2.*synapticChange.D*1000.,model5Hz[:,18]*100.)
        # else:
        #ax2.plot(modelNew[:, 0] * 1000., modelNew[:, 3] * 100.)

        # removes upper and right axes
        # and moves left and bottom axes away
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_position(('outward', 10))
        ax2.spines['left'].set_position(('outward', 10))
        ax2.yaxis.set_ticks_position('left')
        ax2.xaxis.set_ticks_position('bottom')

        ax2.set_xlim(-300, 300)
        # legends and labels
        # plt.legend(loc=1,frameon=False)

        #plt.xlabel(r'$\Delta t$ (ms)')
        plt.ylabel('change in synaptic strength')

        # fourth sub-plot #######################################################
        ax3 = plt.subplot(gs[3])

        # title
        ax3.set_title('10 Hz')

        # diplay of data
        ax3.axhline(y=100, ls='--', color='0.7', lw=2)
        ax3.axvline(x=0, ls='--', color='0.7', lw=2)
        ax3.plot(stdp10Hz[:, 0], stdp10Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5', clip_on=False)
        ax3.errorbar(stdp10binned[:, 0], stdp10binned[:, 1], yerr=stdp10binned[:, 2], fmt='o-', markeredgecolor='C0')
        ax3.plot(synChange[:, 0] * 1000., synChange[:, 4] * 100.)
        if irrData:
            ax3.plot(modelNew[:,0]*1000.,modelNew[:,13]*100./0.5)
            ax3.plot(modelNewReg[:, 0] * 1000., modelNewReg[:, 13] * 100. / 0.5,ls='--')
            # if 'old'==oldNew:
        # ax3.plot(-model10Hz[:,3]+2.*synapticChange.D*1000.,model10Hz[:,18]*100.)
        # else:
        #ax3.plot(modelNew[:, 0] * 1000., modelNew[:, 4] * 100.)

        # removes upper and right axes
        # and moves left and bottom axes away
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['bottom'].set_position(('outward', 10))
        ax3.spines['left'].set_position(('outward', 10))
        ax3.yaxis.set_ticks_position('left')
        ax3.xaxis.set_ticks_position('bottom')

        ax3.set_xlim(-300, 300)
        # legends and labels
        # plt.legend(loc=1,frameon=False)
        if len(paraOpt[0])==7:
            ax3.text(100,230,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s' % (paraOpt[0][0],paraOpt[0][1],paraOpt[0][2],sChange.thetaD,sChange.thetaP,paraOpt[0][3],paraOpt[0][4],paraOpt[0][5],paraOpt[0][6]),fontsize=9)
        elif len(paraOpt[0])==8:
            ax3.text(100,230,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s' % (paraOpt[0][0],paraOpt[0][1],paraOpt[0][2],sChange.thetaD,paraOpt[0][3],paraOpt[0][4],paraOpt[0][5],paraOpt[0][6],paraOpt[0][7]),fontsize=9)

        if not irrData :
            plt.xlabel(r'$\Delta t$ (ms)')
        # plt.ylabel('change in synaptic strength')

        # third sub-plot #######################################################
        ax4 = plt.subplot(gs[4])

        # title
        ax4.set_title('regular : all frequency solutions')

        # diplay of data
        ax4.axhline(y=100, ls='--', color='0.7', lw=2)
        ax4.axvline(x=0, ls='--', color='0.7', lw=2)
        #ax4.plot(stdp5Hz[:, 0], stdp5Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
        #ax4.errorbar(stdp5binned[:, 0], stdp5binned[:, 1], yerr=stdp5binned[:, 2], fmt='o-', markeredgecolor='C0')
        for i in range(len(self.stimFrequencies)):
            ax4.plot(synChange[:, 0] * 1000., synChange[:, i+1] * 100.,label=str(self.stimFrequencies[i]))
        #ax4.plot(synChange[:, 0] * 1000., synChange[:, 2] * 100.)
        #ax4.plot(synChange[:, 0] * 1000., synChange[:, 3] * 100.)

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

        ax4.set_xlim(-300, 300)
        # legends and labels
        plt.legend(loc=1,frameon=False)

        plt.xlabel(r'$\Delta t$ (ms)')
        plt.ylabel('change in synaptic strength')

        # third sub-plot #######################################################
        if irrData :
            ax5 = plt.subplot(gs[5])

            # title
            ax5.set_title('irregular : all frequency solutions')

            # diplay of data
            ax5.axhline(y=100, ls='--', color='0.7', lw=2)
            ax5.axvline(x=0, ls='--', color='0.7', lw=2)
            # ax4.plot(stdp5Hz[:, 0], stdp5Hz[:, 1], 'o', ms=4, c='0.5', markeredgecolor='0.5')
            # ax4.errorbar(stdp5binned[:, 0], stdp5binned[:, 1], yerr=stdp5binned[:, 2], fmt='o-', markeredgecolor='C0')
            for i in range(len(self.stimFrequencies)):
                ax5.plot(modelNew[:,0]*1000.,modelNew[:,(10+i)]*100./0.5)
                #ax5.plot(synChange[:, 0] * 1000., synChange[:, i + 1] * 100., label=str(self.stimFrequencies[i]))
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

            ax5.set_xlim(-300, 300)
            # legends and labels
            plt.legend(loc=1, frameon=False)

            plt.xlabel(r'$\Delta t$ (ms)')
            #plt.ylabel('change in synaptic strength')

        ## save figure ############################################################
        #fname = 'regularDataFitVenance' #_dataSet#' + str(dataSetNumber)

        if irrData:
            fname = 'regular-irregular-DataFitVenance_%s' % paraName #os.path.basename(__file__)
        else:
            fname = 'regularDataFitVenance_%s' % paraName #os.path.basename(__file__)

        savefig(self.figDir + fname + '.png')
        savefig(self.figDir + fname + '.pdf')
        #clf()
        #dsN += 1