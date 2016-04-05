from pylab import *
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import time

from timeAboveThreshold import *

class synUtils(): # synUtils(thetaD,thetaP,nonlinear,Npresentations,w0)
    ###############################################################################
    def __init__(self, thetaD, thetaP, nonlinear, Npresentations, w0):
        self.thetaD  = thetaD
        self.thetaP  = thetaP
        # determine eta based on nonlinearity factor and amplitudes
        self.nonlinear = nonlinear
        self.Npresentations = Npresentations
        self.w0 = w0
        
        # read in experimental data
        dataDir = 'experimental_data/'
        
        jesperReg = loadtxt(dataDir+'sjoestroem_regular_all.dat')
        jesperStoch = loadtxt(dataDir+'sjoestroem_stochastic.dat')

        self.xDataReg = jesperReg[:,[0,1]]
        self.xDataReg[:,1] = self.xDataReg[:,1]/1000. # everything in sec
        self.yDataReg = jesperReg[:,2]+1. # Sjoestroem's data is normalized to 0
        self.sigmaDataReg = jesperReg[:,3]
        
        
        self.xDataStoch = jesperStoch[:,0]
        self.yDataStoch = jesperStoch[:,1]+1. # Sjoestroem's data is normalized to 0
        self.sigmaDataStoch = jesperStoch[:,2]
        
    ##########################################################################################
    # calculate change for regular spike-pair vs frequency protocol
    def calculateChangeInSynapticStrength(self, frequency,deltaT,params):
        #####
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
        tat = timeAboveThreshold(self.thetaD, self.thetaP, tauCa, Cpre, Cpost, self.nonlinear)
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
        mean   =  rhoBar - (rhoBar- 0.5)*exp(-self.Npresentations*interval/tauEff)
        # change in synaptic strength after/before
        return (mean/self.w0)
    
    #############################################################################################
    # calculate change for regular spike-pair vs frequency protocol
    def calculateChangeInSynapticStrengthStochastic(self, frequency,params,DeltaTRange):
        #####
        tauCa = params[0]
        Cpre = params[1]
        Cpost = params[2]
        #thetaD = params[3]
        #thetaP = params[4]
        gammaD = params[3]
        gammaP = params[4]
        tau = params[5]
        D = params[6]
        
        DeltaTStart = DeltaTRange[0]
        DeltaTEnd = DeltaTRange[1]
        
        deltaTs = linspace(DeltaTStart,DeltaTEnd,101)
        
        interval    = 1./frequency
        ####
        tat = timeAboveThreshold(self.thetaD, self.thetaP, tauCa, Cpre, Cpost, self.nonlinear)
        timeDAvg = 0.
        timePAvg = 0.
        for i in range(len(deltaTs)):
            (timeD,timeP) = tat.spikePairFrequencyNonlinear(deltaTs[i]-D,frequency)
            timeDAvg += timeD/len(deltaTs)
            timePAvg += timeP/len(deltaTs)
        GammaP = gammaP*timePAvg/interval
        GammaD = gammaD*timeDAvg/interval
        # rhoBar: average value of rho in the limit of a very long protocol equivalent to the minimum of the quadratic potentia
        try :
            rhoBar = GammaP/(GammaP + GammaD)
        except RuntimeWarning:
            print GammaP, GammaD
        # tauEff : characteristic time scale of the temporal evolution of the pdf of rho
        tauEff = tau/(GammaP + GammaD)
        #
        # mean value of the synaptic strength right at the end of the stimulation protocol
        mean   =  rhoBar - (rhoBar- 0.5)*exp(-self.Npresentations*interval/tauEff)
        # change in synaptic strength after/before
        return (mean/self.w0)
    
    ################################################################################################
    def generateFig(self, paraOpt):
        
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
        ax0.set_title('regular Sjoestroem , chi2='+str(paraOpt[1]))

        # diplay of data
        ax0.axhline(y=1.,c='0.7')
        ax0.plot(self.xDataReg[:,0][::2],self.yDataReg[::2],'s',color='red',clip_on=False)
        ax0.plot(self.xDataReg[:,0][1::2],self.yDataReg[1::2],'o',color='blue',clip_on=False)
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
        ax0.plot(self.xDataStoch,self.yDataStoch,'s',color='green',clip_on=False)
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
        fname = 'synChangeCa-ModelFit' #os.path.basename(__file__)

        savefig(fname[:-3]+'.png')
        savefig(fname[:-3]+'.pdf')
