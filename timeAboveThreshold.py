#import scipy as sci
import numpy as np
import pdb
import sys
#import commands
#import random as rrr

class timeAboveThreshold():
        ''' 
            class to calculate the fraction of time \alpha the calcium trace spends above threshold 
        '''
        ###############################################################################
        def __init__(self, tauCa, Cpre, Cpost, thetaD, thetaP, nonlinear=1.,Nves=0):
                self.tauCa = tauCa
                self.Cpre = Cpre
                self.Cpost = Cpost
                self.thetaD  = thetaD
                self.thetaP  = thetaP
                # determine eta based on nonlinearity factor and amplitudes
                self.nonlinear = nonlinear
                self.eta = (self.nonlinear*(self.Cpost + self.Cpre) - self.Cpost)/self.Cpre - 1.
                self.Nvesicles = Nves
                
        ###############################################################################
        # regular spike-pairs vs. frequency
        def spikePairFrequency(self,deltaT,frequency):
                #
                
                interval = 1./frequency
                
                timeAbove = zeros(2)

                # in case deltaT is larger then one interval
                while (deltaT > 1. / (2. * frequency)):
                        deltaT = deltaT - 1. / frequency
                while (deltaT < -1. / (2. * frequency)):
                        deltaT = deltaT + 1. / frequency

                # determine amplitude of the discontinous points of the calcium trace
                # post-pre
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        A = 0.
                        B = self.Cpost*exp(-fabs(deltaT)/self.tauCa)
                else :
                        A  = (self.Cpost + self.Cpre*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        B  = (self.Cpre + self.Cpost*exp((1./frequency -fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                C  = A + self.Cpost
                D  = B + self.Cpre

                # pre-post
                if ( exp(1./(frequency*self.tauCa)) == NaN ) :
                        E = 0.
                        F = self.Cpre*exp(-fabs(deltaT)/self.tauCa)
                else :
                        E  = (self.Cpre + self.Cpost*exp(fabs(deltaT)/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        F  = (self.Cpost + self.Cpre*exp((1./frequency -fabs(deltaT))/self.tauCa))/(exp(1./(frequency*self.tauCa)) - 1.)
                        
                G  = E + self.Cpre
                H  = F + self.Cpost
                
                # loop over depression and potentiation threshold
                for i in range(2):
                        if i == 0:
                                Ct = self.thetaD
                        elif i==1:
                                Ct = self.thetaP
                        # post-pre
                        if (deltaT < 0.):
                                if ( A <= Ct and B > Ct ) :
                                        I = self.tauCa*np.log(D/Ct) + fabs(deltaT)
                                elif ( A > Ct and B > Ct) :
                                        I = 1./frequency
                                elif ( A > Ct and B <= Ct ) :
                                        I = self.tauCa*np.log(C/Ct) + fabs(1./frequency) - fabs(deltaT)
                                elif ( A <= Ct and B <= Ct and D > Ct and C > Ct ) :
                                        I = self.tauCa*np.log(C/Ct) + self.tauCa*np.log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C > Ct ) :
                                        I = self.tauCa*np.log(C/Ct)
                                elif( A <= Ct and B <= Ct and D > Ct and C <= Ct ) :
                                        I = self.tauCa*np.log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C <= Ct ) :
                                        I = 0.
                                else :
                                        print(A, B, C, D, Ct, frequency, deltaT)
                                        print("post-pre : Problem in spikePairFrequency!")
                                        sys.exit(1)
                        # pre-post
                        else:
                                if ( E <= Ct and F > Ct ) :
                                        I = self.tauCa*np.log(H/Ct) + fabs(deltaT)
                                elif ( E > Ct and F > Ct) :
                                        I = 1./frequency
                                elif ( E > Ct and F <= Ct ) :
                                        I = self.tauCa*np.log(G/Ct) + fabs(1./frequency) - fabs(deltaT)
                                elif ( E <= Ct and F <= Ct and G > Ct and H > Ct ) :
                                        I = self.tauCa*np.log(G/Ct) + self.tauCa*np.log(H/Ct)
                                elif ( E <= Ct and F <= Ct and H <= Ct and G > Ct ) :
                                        I = self.tauCa*np.log(G/Ct)
                                elif ( E <= Ct and F <= Ct and H > Ct and G <= Ct ) :
                                        I = self.tauCa*np.log(H/Ct)
                                elif ( E <= Ct and F <= Ct and G <= Ct and H <= Ct ) :
                                        I = 0.
                                else :
                                        print(E, F, G, H, Ct, frequency, deltaT)
                                        print("pre-post : Problem in spikePairFrequency! ")
                                        sys.exit(1)
                        #
                        timeAbove[i] = I
                #
                alphaD = timeAbove[0]/interval
                alphaP = timeAbove[1]/interval
                return (alphaD,alphaP)
        ###############################################################################
        # regular spike-pairs at a given frequency, implements the nonlinear calcium model
        def spikePairFrequencyNonlinear(self,deltaT,frequency):
                #
                interval = 1./frequency
                
                timeAbove = np.zeros(2)
                
                # in case deltaT is larger then one interval
                while (deltaT > 1./(2.*frequency) ):
                        deltaT = deltaT - 1./frequency
                while (deltaT < -1./(2.*frequency) ):
                        deltaT = deltaT + 1./frequency


                # determine amplitude of the discontinous points of the calcium trace
                # post-pre
                if ( np.exp(1./(frequency*self.tauCa)) == np.NaN ) :
                        A = 0.
                        B = self.Cpost*np.exp(-fabs(deltaT)/self.tauCa)
                else :
                        A  = (self.Cpost + (self.eta*np.exp(-1./(frequency*self.tauCa))/(1.-np.exp(-1./(frequency*self.tauCa))) + 1.)*self.Cpre*np.exp(np.fabs(deltaT)/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) - 1.)
                        B  = self.Cpre*(1. + self.eta/(1.-np.exp(-1./(frequency*self.tauCa))))/(np.exp(1./(frequency*self.tauCa)) - 1.) + self.Cpost*np.exp(-np.fabs(deltaT)/self.tauCa)/(1. - np.exp(-1./(frequency*self.tauCa)))
                C  = A + self.Cpost + self.eta*self.Cpre*np.exp((np.fabs(deltaT)-1./frequency)/self.tauCa)/(1.-np.exp(-1./(frequency*self.tauCa)))
                D  = B + self.Cpre

                # pre-post
                if ( np.exp(1./(frequency*self.tauCa)) == np.NaN ) :
                        E = 0.
                        F = self.Cpre*np.exp(-np.fabs(deltaT)/self.tauCa)
                else :
                        E  = (self.Cpre*(1. + self.eta/(1.-np.exp(-1./(frequency*self.tauCa)))) + self.Cpost*np.exp(np.fabs(deltaT)/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) - 1.)
                        F  = (self.Cpost + self.eta/(1.-np.exp(-1./(frequency*self.tauCa)))*self.Cpre*np.exp(-np.fabs(deltaT)/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) - 1.) + self.Cpre*np.exp(-np.fabs(deltaT)/self.tauCa)/(1. - np.exp(-1./(frequency*self.tauCa)))
                        
                G  = E + self.Cpre 
                H  = F + self.Cpost + self.eta*self.Cpre*np.exp(-np.fabs(deltaT)/self.tauCa)/(1.-np.exp(-1./(frequency*self.tauCa)))
                
                # loop over depression and potentiation threshold
                for i in range(2):
                        if i == 0:
                                Ct = self.thetaD
                        elif i==1:
                                Ct = self.thetaP
                        # post-pre
                        if (deltaT < 0.):
                                if ( A <= Ct and B > Ct ) :
                                        I = self.tauCa*np.log(D/Ct) + np.fabs(deltaT)
                                elif ( A > Ct and B > Ct) :
                                        I = 1./frequency
                                elif ( A > Ct and B <= Ct ) :
                                        I = self.tauCa*np.log(C/Ct) + np.fabs(1./frequency) - np.fabs(deltaT)
                                elif ( A <= Ct and B <= Ct and D > Ct and C > Ct ) :
                                        I = self.tauCa*np.log(C/Ct) + self.tauCa*np.log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C > Ct ) :
                                        I = self.tauCa*np.log(C/Ct)
                                elif( A <= Ct and B <= Ct and D > Ct and C <= Ct ) :
                                        I = self.tauCa*np.log(D/Ct)
                                elif ( A <= Ct and B <= Ct and D <= Ct and C <= Ct ) :
                                        I = 0.
                                else :
                                        print(A, B, C, D, Ct, frequency, deltaT)
                                        print("post-pre : Problem in spikePairFrequency!")
                                        sys.exit(1)
                        # pre-post
                        else:
                                if ( E <= Ct and F > Ct ) :
                                        I = self.tauCa*np.log(H/Ct) + np.fabs(deltaT)
                                elif ( E > Ct and F > Ct) :
                                        I = 1./frequency
                                elif ( E > Ct and F <= Ct ) :
                                        I = self.tauCa*np.log(G/Ct) + np.fabs(1./frequency) - np.fabs(deltaT)
                                elif ( E <= Ct and F <= Ct and G > Ct and H > Ct ) :
                                        I = self.tauCa*np.log(G/Ct) + self.tauCa*np.log(H/Ct)
                                elif ( E <= Ct and F <= Ct and H <= Ct and G > Ct ) :
                                        I = self.tauCa*np.log(G/Ct)
                                elif ( E <= Ct and F <= Ct and H > Ct and G <= Ct ) :
                                        I = self.tauCa*np.log(H/Ct)
                                elif ( E <= Ct and F <= Ct and G <= Ct and H <= Ct ) :
                                        I = 0.
                                else :
                                        print(E, F, G, H, Ct, frequency, deltaT)
                                        print("pre-post : Problem in spikePairFrequency! ")
                                        sys.exit(1)
                        #
                        timeAbove[i] = I
                #
                alphaD = timeAbove[0]/interval
                alphaP = timeAbove[1]/interval
                return (alphaD,alphaP)
        
        ###############################################################################
        # one presynaptic spike and a two postsynaptic spikes, i.e., a postsynaptic burst
        def preSpikePostPair(self,deltaT,frequency,deltaBurst):
                #
                interval = 1./frequency
                
                timeAbove = zeros(2)
                #
                if ( np.fabs(deltaT) > 1./frequency ):
                        deltaT = -(np.fabs(deltaT) - 1./frequency)
                        
                #########################################
                # post-post-pre
                M = self.Cpost*np.exp(-(np.fabs(deltaBurst))/self.tauCa)
                N = self.Cpost*np.exp(-(np.fabs(deltaBurst) + np.fabs(deltaT))/self.tauCa) + self.Cpost*np.exp(-np.fabs(deltaT)/self.tauCa)

                if ( np.exp(1./(frequency*self.tauCa)) == NaN ):
                        A = 0.
                else :
                        A  = (self.Cpost + self.Cpost*np.exp(np.fabs(deltaBurst)/self.tauCa) + self.Cpre*np.exp((np.fabs(deltaBurst) + np.fabs(deltaT))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)
                        M = (self.Cpost + self.Cpre*np.exp(np.fabs(deltaT)/self.tauCa) + self.Cpost*np.exp((1./frequency - np.fabs(deltaBurst))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)
                        N = (self.Cpre + self.Cpost*np.exp((1./frequency - np.fabs(deltaT) - np.fabs(deltaBurst))/self.tauCa) + self.Cpost*np.exp((1./frequency - np.fabs(deltaT))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)

                O  = M + self.Cpost
                P  = N + self.Cpre
                
                #########################################
                # post-pre-post
                Q = self.Cpost*np.exp(-(np.fabs(deltaBurst) - np.fabs(deltaT))/self.tauCa)
                R = self.Cpost*np.exp(- np.fabs(deltaBurst)/self.tauCa) + self.Cpre*np.exp(-np.fabs(deltaT)/self.tauCa)
                if ( np.exp(1./(frequency*self.tauCa)) == NaN ) :
                        A = 0.
                else :
                        A  = (self.Cpost + self.Cpre*np.exp((np.fabs(deltaBurst)-np.fabs(deltaT))/self.tauCa) + self.Cpost*np.exp(np.fabs(deltaBurst)/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)
                        Q = (self.Cpre + self.Cpost*np.exp(np.fabs(deltaT)/self.tauCa) + self.Cpost*np.exp((1./frequency - (np.fabs(deltaBurst) - np.fabs(deltaT)))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)
                        R = (self.Cpost + self.Cpost*np.exp((1./frequency - np.fabs(deltaBurst))/self.tauCa) + self.Cpre*np.exp((1./frequency - np.fabs(deltaT))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)

                S  = Q + self.Cpre
                T  = R + self.Cpost
                
                
                #########################################
                # pre-post-post
                U = self.Cpre*np.exp(-(np.fabs(deltaT)-np.fabs(deltaBurst))/self.tauCa)
                V = self.Cpre*np.exp(- np.fabs(deltaT)/self.tauCa) + self.Cpost*np.exp(-np.fabs(deltaBurst)/self.tauCa)
                
                if ( np.exp(1./(frequency*self.tauCa)) == NaN ):
                        A = 0.
                else :
                        A  = (self.Cpre + self.Cpost*np.exp((np.fabs(deltaT)-np.fabs(deltaBurst))/self.tauCa) + self.Cpost*np.exp(np.fabs(deltaT)/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)
                        U = (self.Cpost + self.Cpost*np.exp(np.fabs(deltaBurst)/self.tauCa) + self.Cpre*np.exp((1./frequency - (np.fabs(deltaT) - np.fabs(deltaBurst)))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)
                        V = (self.Cpost + self.Cpre*np.exp((1./frequency - np.fabs(deltaT))/self.tauCa) + self.Cpost*np.exp((1./frequency - np.fabs(deltaBurst))/self.tauCa))/(np.exp(1./(frequency*self.tauCa)) -1.)

                W  = U + self.Cpost
                X  = V + self.Cpost
	

                for i in range(2):
                        if i == 0:
                                Ct = self.thetaD
                        elif i==1:
                                Ct = self.thetaP
                        # 2 x post and 1 x pre	
                        # post-post-pre
                        if ( deltaT < 0.) :
                                # print "post-post-pre", deltaT, deltaBurst
                                if ( M > Ct and N > Ct and A > Ct) :
                                        Int = 1./frequency
                                elif ( M > Ct and N > Ct and A <= Ct) :
                                        Int = self.tauCa*np.log(P/Ct) + np.fabs(deltaT) + np.fabs(deltaBurst)
                                elif ( M > Ct and N <= Ct and P > Ct) :
                                        Int = self.tauCa*np.log(O/Ct) + np.fabs(deltaBurst) + self.tauCa*np.log(P/Ct)
                                elif ( M > Ct and P <= Ct) :
                                        Int = self.tauCa*np.log(O/Ct) + np.fabs(deltaBurst)
                                elif ( (A+self.Cpost) > Ct and M <= Ct and N > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(P/Ct) + np.fabs(deltaT)
                                elif ( (A+self.Cpost) > Ct and M <= Ct and N <= Ct and P > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(O/Ct) +  self.tauCa*np.log(P/Ct)
                                elif ( (A+self.Cpost) > Ct and M <= Ct and P <= Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(O/Ct)
                                elif ( (A+self.Cpost) > Ct and O <= Ct and P <= Ct) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct)
                                elif ( (A+self.Cpost) > Ct and O <= Ct and P > Ct) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(P/Ct)
                                elif ( (A+self.Cpost) <= Ct and N > Ct ) :
                                        Int = self.tauCa*np.log(P/Ct) + np.fabs(deltaT)
                                elif ( (A+self.Cpost) <= Ct and O > Ct and N <= Ct and P > Ct) :
                                        Int = self.tauCa*np.log(O/Ct) + self.tauCa*np.log(P/Ct)
                                elif ( (A+self.Cpost) <= Ct and O > Ct and P <= Ct) :
                                        Int = self.tauCa*np.log(O/Ct)
                                elif ( (A+self.Cpost) <= Ct and O <= Ct and P > Ct) :
                                        Int = self.tauCa*np.log(P/Ct)
                                elif ( (A+self.Cpost) <= Ct and O <= Ct and P <= Ct) :
                                        Int = 0.
                                else :
                                        print("post-post-pre : Problem !")
                                        print(deltaT, A, (A+self.Cpost), M, N, O, P, Int)
                                        sys.exit(1)
                        # post-pre-post 
                        elif ( deltaT >= 0.  and deltaT <= deltaBurst )  :
                                # print "post-pre-post", deltaT, deltaBurst 
                                if ( Q > Ct and R > Ct and A > Ct) :
                                        Int = 1./frequency
                                elif ( Q > Ct and R > Ct and A <= Ct) :
                                        Int = self.tauCa*np.log(T/Ct) + np.fabs(deltaBurst)
                                elif ( Q > Ct and R <= Ct and T > Ct ) :
                                        Int = self.tauCa*np.log(S/Ct) + np.fabs(deltaBurst) - np.fabs(deltaT) +  self.tauCa*np.log(T/Ct)
                                elif ( Q > Ct and T <= Ct and A <= Ct ) :
                                        Int = self.tauCa*np.log(S/Ct) + np.fabs(deltaBurst) - np.fabs(deltaT)
                                elif ( (A+self.Cpost) > Ct and Q <= Ct and R > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(T/Ct) + np.fabs(deltaT)
                                elif ( (A+self.Cpost) > Ct and Q <= Ct and S > Ct and R <= Ct and T > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(S/Ct) + self.tauCa*np.log(T/Ct)
                                elif ( (A+self.Cpost) > Ct and Q <= Ct and S > Ct and T <= Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(S/Ct)
                                elif ( (A+self.Cpost) > Ct and S <= Ct and T > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct) + self.tauCa*np.log(T/Ct)
                                elif ( (A+self.Cpost) > Ct and S <= Ct and T <= Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpost)/Ct)
                                elif ( (A+self.Cpost) <= Ct and S > Ct and T <= Ct ) :
                                        Int = self.tauCa*np.log(S/Ct)
                                elif ( (A+self.Cpost) <= Ct and S <= Ct and T > Ct) :
                                        Int = self.tauCa*np.log(T/Ct)
                                elif ( (A+self.Cpost) <= Ct and S > Ct and R <= Ct and T > Ct) :
                                        Int = self.tauCa*np.log(S/Ct) + self.tauCa*np.log(T/Ct)
                                elif ( (A+self.Cpost) <= Ct  and R > Ct) :
                                        Int = self.tauCa*np.log(T/Ct) + np.fabs(deltaT)
                                elif ( (A+self.Cpost) <= Ct and S <= Ct and T <= Ct) :
                                        Int = 0.
                                else :
                                        print("post-pre-post : Problem !")
                                        print(deltaT, A, Q, R, (A+self.Cpost), S, T, Int)
                                        sys.exit(1)
                        # pre-post-post 
                        elif ( deltaT >= 0.  and deltaT > deltaBurst)  :
                                # print "pre-post-post", deltaT, deltaBurst
                                if ( U > Ct and V > Ct and A > Ct ) :
                                        Int = 1./frequency
                                elif ( U > Ct and V > Ct and A <= Ct) :
                                        Int = self.tauCa*np.log(X/Ct) + np.fabs(deltaT)
                                elif ( U > Ct and V <= Ct and X > Ct ) :
                                        Int = self.tauCa*np.log(W/Ct) + np.fabs(deltaT) - np.fabs(deltaBurst) + self.tauCa*np.log(X/Ct)
                                elif ( U > Ct and X <= Ct ) :
                                        Int = self.tauCa*np.log(W/Ct) + np.fabs(deltaT) - np.fabs(deltaBurst)
                                elif ( (A+self.Cpre) > Ct and U <= Ct and  V > Ct) :
                                        Int = self.tauCa*np.log((A+self.Cpre)/Ct) + self.tauCa*np.log(X/Ct) + np.fabs(deltaBurst)
                                elif ( (A+self.Cpre) > Ct and U <= Ct and  V <= Ct and W > Ct and X > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpre)/Ct) + self.tauCa*np.log(W/Ct) + self.tauCa*np.log(X/Ct)
                                elif ( (A+self.Cpre) > Ct and U <= Ct and W > Ct and X <= Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpre)/Ct) + self.tauCa*np.log(W/Ct)
                                elif ( (A+self.Cpre) > Ct and W <= Ct and X > Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpre)/Ct) + self.tauCa*np.log(X/Ct)
                                elif ( (A+self.Cpre) > Ct and W <= Ct and X <= Ct ) :
                                        Int = self.tauCa*np.log((A+self.Cpre)/Ct)
                                elif ( (A+self.Cpre) <= Ct and V > Ct ) :
                                        Int = self.tauCa*np.log(X/Ct) + np.fabs(deltaBurst)
                                elif ( (A+self.Cpre) <= Ct  and V <= Ct and W > Ct and X > Ct ) :
                                        Int = self.tauCa*np.log(W/Ct) + self.tauCa*np.log(X/Ct)
                                elif ( (A+self.Cpre) <= Ct  and W <= Ct and X > Ct ) :
                                        Int = self.tauCa*np.log(X/Ct)
                                elif ( (A+self.Cpre) <= Ct  and W > Ct and X <= Ct ) :
                                        Int = self.tauCa*np.log(W/Ct)
                                elif ( (A+self.Cpre) <= Ct  and W <= Ct and X <= Ct ) :
                                        Int = 0.
                                else :
                                        print("pre-post-post : Problem !")
                                        print(deltaT, A, (A+self.Cpre), U, V, W, X, Int)
                                        sys.exit(1)
                        else :
                                print("Problem in preSpikePostPair routine!")
                                sys.exit(1)
                        #
                        timeAbove[i] = Int
                #
                alphaD = timeAbove[0]/interval
                alphaP = timeAbove[1]/interval
                return (alphaD,alphaP)
        ###############################################################################
        # irregular spike-pairs, the numerical integration is run in an external C++ code for performance improvment
        def irregularSpikePairs(self,deltaT,preRate,postRate,ppp,deltaCa):
                
                # linear calcium dynamics, numerical integration possible
                if self.nonlinear == 1.:
                    #print 'time above threshold : integration for LINEAR calcium dynamics'
                    # the first argument calcium ampliutde has to be smaller than the second
                    if self.Cpre>self.Cpost:
                            arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpost) + ' ' + str(self.Cpre) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
                    else:
                            arguments = str(deltaT) + ' ' + str(self.tauCa) + ' ' + str(self.Cpre) + ' ' + str(self.Cpost) + ' ' + str(self.thetaD) + ' ' + str(self.thetaP) + ' ' + str(preRate) + ' ' + str(postRate) + ' ' + str(ppp) + ' ' + str(deltaCa)
        
                    #print arguments
                    #(out,err) = commands.getstatusoutput('./timeAboveThreshold/poissonPairs_timeAboveThreshold ' + arguments)
                    alphaD = float(err.split()[0])
                    alphaP = float(err.split()[1])
                    
                # nonlinear calcium dynamics
                else:
                    #print 'time above threshold : NONLINEAR calcium dynamics'
                    
                    # construction of the spike train
                    #tStart = 0.1 # start time at 100 ms
                    
                    rrr.seed(7)
                    
                    tPre = []
                    tPostCorr = []
                    tPostInd = []
                    
                    tPre.append(0)
                    tPostInd.append(0)
                    
                    for i in range(10000):
                        tPre.append(tPre[-1] + np.random.exponential(1./preRate))
                        if np.random.rand()<ppp:
                            tPostCorr.append(tPre[-1]+deltaT)
                        if (postRate-ppp*preRate) > 0.:
                            tPostInd.append(tPostInd[-1] + np.random.exponential(1./(postRate-ppp*preRate)))
                    
                    
                    tPost = tPostCorr + tPostInd[1:]
                    
                    tPostSorted = sorted(tPost, key=lambda tPost: tPost)
                    
                    tAll = zeros((len(tPre[1:]) + len(tPostSorted),3))
                    
                    tAll[:,0] = hstack((tPre[1:],tPostSorted))
                    tAll[:,1] = hstack((zeros(len(tPre[1:])),ones(len(tPostSorted))))
                    tAll[:,2] = hstack((repeat(self.Cpre,len(tPre[1:]),repeat(self.Cpost,len(tPostSorted)))))

                    tList = tAll.tolist()
                    tListSorted = sorted(tList, key=lambda tList: tList[0])
                    
                    #tListSorted.append([Npres/freq,2])

                    (tD, tP) = self.eventBasedIntegration(tListSorted)

                    alphaD = tD/tListSorted[-1][0]
                    alphaP = tP/tListSorted[-1][0]
                    #print alphaD, alphaP
                ####################################################################
                
                return (alphaD,alphaP)

        ###############################################################################
        # irregular spike-pairs, the numerical integration is run in an external C++ code for performance improvment
        def irregularSpikePairsEventBased(self, deltaT, preRate, postRate, ppp, nSpikes=10000):


            # print 'time above threshold : NONLINEAR calcium dynamics'

            # construction of the spike train
            # tStart = 0.1 # start time at 100 ms

            #np.random.seed(7)

            tPre = []
            tPostCorr = []
            tPostInd = []

            tPre.append(0)
            tPostInd.append(0)

            for i in range(nSpikes):
                tPre.append(tPre[-1] + np.random.exponential(1. / preRate))
                if np.random.rand() < ppp:
                    tPostCorr.append(tPre[-1] + deltaT)
                if (postRate - ppp * preRate) > 0.:
                    tPostInd.append(tPostInd[-1] + np.random.exponential(1. / (postRate - ppp * preRate)))

            tPost = tPostCorr + tPostInd[1:]

            tPostSorted = sorted(tPost, key=lambda tPost: tPost)

            tAll = np.zeros((len(tPre[1:]) + len(tPostSorted), 3))

            tAll[:, 0] = np.hstack((tPre[1:], tPostSorted))
            tAll[:, 1] = np.hstack((np.zeros(len(tPre[1:])), np.ones(len(tPostSorted))))
            tAll[:, 2] = np.hstack((np.repeat(self.Cpre, len(tPre[1:])), np.repeat(self.Cpost, len(tPostSorted))))

            tList = tAll.tolist()
            tListSorted = sorted(tList, key=lambda tList: tList[0])

            # tListSorted.append([Npres/freq,2])

            (tD, tP) = self.eventBasedIntegration(tListSorted)

            alphaD = tD / tListSorted[-1][0]
            alphaP = tP / tListSorted[-1][0]  # print alphaD, alphaP
            ####################################################################

            return (alphaD, alphaP)

        ############################################################################################
        ## single out bursts from plasticity trace #################################################
        def separateBursts(self,spikeTimes,burstInterval):
            onlyBurstsTemp = []
            onlyBursts = []
            noBursts =  []
            #onlyBursts.append(spikeTimes[0])
            noBursts.append(spikeTimes[0])
            for i in range(1,len(spikeTimes)):
                if (spikeTimes[i] - spikeTimes[i-1])<burstInterval:
                    onlyBurstsTemp.append(spikeTimes[i-1])
                    onlyBurstsTemp.append(spikeTimes[i])
                else:
                    noBursts.append(spikeTimes[i])


            # remove duplicates
            for i in onlyBurstsTemp:
                if i not in onlyBursts:
                    onlyBursts.append(i)

            onlyBursts = np.asarray(onlyBursts)
            noBursts = np.asarray(noBursts)
            # count bursts and number of spikes in bursts
            diffB = onlyBursts[1:]-onlyBursts[:-1]
            numberOfBursts = np.sum(diffB>burstInterval) + 1 # number of bursts is given by intervals larger than the 'burstInterval'
            boolBursts = diffB<burstInterval
            boolBursts = np.concatenate((np.array([False]),boolBursts))
            boolBursts = np.concatenate((np.diff(boolBursts),np.array([True])))
            startStopBursts = np.arange(len(boolBursts))[boolBursts]
            #[i for i,(m,n) in enumerate(zip([2]+boolBursts,boolBursts+[2])) if m!=n]
            #startStopBursts = np.asarray(startStopBursts)
            numberSpikesInBursts = (startStopBursts[1::2]-startStopBursts[:-1:2])+1
            #pdb.set_trace()
            return (onlyBursts,noBursts,numberOfBursts,numberSpikesInBursts)

        ###############################################################################
        # irregular spike-pairs, bursts only, the numerical integration is run in an external C++ code for performance improvment
        def irregularBurstPairsEventBased(self, deltaT, preRate, postRate, ppp,nSpikes=2000):

                # print 'time above threshold : NONLINEAR calcium dynamics'

                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                #np.random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(nSpikes):
                        tPre.append(tPre[-1] + np.random.exponential(1. / preRate))
                        if np.random.rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(tPostInd[-1] + np.random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)

                (tPreBursts,_,numberOfBursts,numberSpikesInBursts)  = self.separateBursts(tPre[1:],0.15)
                (tPostBursts,_,numberOfBursts,numberSpikesInBursts) = self.separateBursts(tPostSorted,0.15)
                #print numberOfBursts, numberSpikesInBursts
                tAll = np.zeros((len(tPreBursts) + len(tPostBursts), 3))

                tAll[:, 0] = np.hstack((tPreBursts, tPostBursts))
                tAll[:, 1] = np.hstack((np.zeros(len(tPreBursts)), np.ones(len(tPostBursts))))
                tAll[:, 2] = np.hstack((np.repeat(self.Cpre, len(tPreBursts)), np.repeat(self.Cpost, len(tPostBursts))))

                tList = tAll.tolist()
                tListSorted = sorted(tList, key=lambda tList: tList[0])

                # tListSorted.append([Npres/freq,2])

                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD / tListSorted[-1][0]
                alphaP = tP / tListSorted[-1][0]  # print alphaD, alphaP
                normalizeF = nSpikes/100.
                return (alphaD, alphaP,numberOfBursts/normalizeF,numberSpikesInBursts)

        ###############################################################################
        # irregular spike-pairs, bursts only, the numerical integration is run in an external C++ code for performance improvment
        def irregularIndSpikePairsEventBased(self, deltaT, preRate, postRate, ppp,nSpikes=2000):

                # print 'time above threshold : NONLINEAR calcium dynamics'

                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                #np.random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(nSpikes):
                        tPre.append(tPre[-1] + np.random.exponential(1. / preRate))
                        if np.random.rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(tPostInd[-1] + np.random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)

                (_,tPreNoBursts, _, _) = self.separateBursts(tPre[1:], 0.15)
                (_, tPostNoBursts, _, _) = self.separateBursts(tPostSorted, 0.15)

                tAll = np.zeros((len(tPreNoBursts) + len(tPostNoBursts), 3))

                tAll[:, 0] = np.hstack((tPreNoBursts, tPostNoBursts))
                tAll[:, 1] = np.hstack((np.zeros(len(tPreNoBursts)), np.ones(len(tPostNoBursts))))
                tAll[:, 2] = np.hstack((np.repeat(self.Cpre, len(tPreNoBursts)), np.repeat(self.Cpost, len(tPostNoBursts))))

                tList = tAll.tolist()
                tListSorted = sorted(tList, key=lambda tList: tList[0])

                # tListSorted.append([Npres/freq,2])

                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD / tListSorted[-1][0]
                alphaP = tP / tListSorted[-1][0]  # print alphaD, alphaP

                return (alphaD, alphaP)
        ###############################################################################
        # irregular spike-pairs and deterministic short-term plasticity, event-based integration
        def irregularSpikePairsSTPDeterministic(self, deltaT, preRate, postRate, ppp, tauRec, U):


                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(50000):
                        tPre.append(tPre[-1] + random.exponential(1. / preRate))
                        if np.random.rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(
                                        tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)

                cpre = zeros(len(tPre[1:]))
                # deterministic STD model implementation
                if U != 0.:
                        cpre[0] = 1.
                        for i in range(1, len(tPre[1:])):
                                cpre[i] = 1. - (1. - (cpre[i - 1] - U * cpre[i - 1])) * np.exp(-(tPre[i+1]-tPre[i])/tauRec)
                        cpre *= U * self.Cpre
                else:
                        cpre[:] = self.Cpre


                tAll = zeros((len(tPre[1:]) + len(tPostSorted), 3))

                tAll[:, 0] = hstack((tPre[1:], tPostSorted))
                tAll[:, 1] = hstack((zeros(len(tPre[1:])), ones(len(tPostSorted))))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,len(tPostSorted))))
                tList = tAll.tolist()
                tListSorted = sorted(tList, key=lambda tList: tList[0])

                # tListSorted.append([Npres/freq,2])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD / tListSorted[-1][0]
                alphaP = tP / tListSorted[-1][0]
                # print alphaD, alphaP

                return (alphaD, alphaP)

        ###############################################################################
        # irregular spike-pairs and deterministic short-term plasticity, event-based integration
        def irregularSpikePairsSTPStochastic(self, deltaT, preRate, postRate, ppp, tauRec, pRelease, Nves):

                NpreSpikes = 5000

                NrepetitionsStoch = 100

                q = self.Cpre/Nves

                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(NpreSpikes):
                        tPre.append(tPre[-1] + random.exponential(1. / preRate))
                        if np.random.rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(
                                        tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)
                tPre = asarray(tPre[1:])

                #######################################
                # stochastic STD model implementation
                #ampStoch = np.zeros((NrepetitionsStoch, len(tPre)))
                timesAbove = zeros((NrepetitionsStoch, 2))
                #tP = 0.
                for r in range(NrepetitionsStoch):
                        Vesicles = np.ones((len(tPre), Nves))
                        Release = np.random.rand(len(tPre), Nves) < pRelease
                        #VesTimes = transpose(np.tile(tPre, (Nves, 1)))
                        for i in range(len(tPre)):
                                nRel = sum(Release[i])
                                releaseSites = argwhere(Release[i] == True)
                                emptyTimes = random.exponential(tauRec, nRel)
                                for n in range(nRel):
                                        if not Vesicles[i, releaseSites[n]] == 0.:
                                                recoveryMask = ((tPre - tPre[i]) < emptyTimes[n]) & ((tPre - tPre[i]) > 0)
                                                Vesicles[recoveryMask, releaseSites[n]] = 0
                        #print 'after amp det'
                        cpreStoch = q*(sum(Vesicles*Release,axis=1))

                        ######################################
                        tAll = zeros((len(tPre) + len(tPostSorted), 3))

                        tAll[:, 0] = hstack((tPre, tPostSorted))
                        tAll[:, 1] = hstack((zeros(len(tPre)), ones(len(tPostSorted))))
                        tAll[:, 2] = hstack((cpreStoch,repeat(self.Cpost,len(tPostSorted))))
                        tList = tAll.tolist()
                        tListSorted = sorted(tList, key=lambda tList: tList[0])

                        # tListSorted.append([Npres/freq,2])

                        ###########################################################
                        # event-based integration
                        (timesAbove[r,0], timesAbove[r,1]) = self.eventBasedIntegration(tListSorted)
                        #print r

                alphaD = average(timesAbove[:,0]) / tListSorted[-1][0]
                alphaP = average(timesAbove[:,1]) / tListSorted[-1][0]
                # print alphaD, alphaP

                return (alphaD, alphaP)


        ###############################################################################
        # irregular spike-pairs and short-term plasticity, event-based integration the numerical integration is run in an external C++ code for performance improvment
        def irregularSpikePairsTMM(self, deltaT, preRate, postRate, ppp, U, tauFac, tauDep):


                # construction of the spike train
                # tStart = 0.1 # start time at 100 ms

                random.seed(7)

                tPre = []
                tPostCorr = []
                tPostInd = []

                tPre.append(0)
                tPostInd.append(0)

                for i in range(10000):
                        tPre.append(tPre[-1] + random.exponential(1. / preRate))
                        if np.random.rand() < ppp:
                                tPostCorr.append(tPre[-1] + deltaT)
                        if (postRate - ppp * preRate) > 0.:
                                tPostInd.append(
                                        tPostInd[-1] + random.exponential(1. / (postRate - ppp * preRate)))

                tPost = tPostCorr + tPostInd[1:]

                tPostSorted = sorted(tPost, key=lambda tPost: tPost)

                cpre = zeros(len(tPre[1:]))
                x    = zeros(len(tPre[1:]))
                u    = zeros(len(tPre[1:]))
                if U != 0.:
                        u[0] = U
                        x[0] = 1.
                        cpre[0] = u[0] * x[0]
                        for i in range(1, len(tPre[1:])):
                                u[i] = u[i - 1] * np.exp(-(tPre[i+1]-tPre[i]) / tauFac) + U * (1. - u[i - 1] * np.exp(-(tPre[i+1]-tPre[i]) / tauFac))
                                x[i] = 1. - (1. - (x[i - 1] - u[i - 1] * x[i - 1])) * np.exp(-(tPre[i+1]-tPre[i]) / tauDep)
                                cpre[i] = u[i] * x[i]
                        cpre *= self.Cpre
                        #cpre[0] = 1.
                        #for i in range(1, len(tPre[1:])):
                        #        cpre[i] = 1. - (1. - (cpre[i - 1] - U * cpre[i - 1])) * np.exp(-(tPre[i+1]-tPre[i])/tauRec)
                        #cpre *= U * self.Cpre
                else:
                        cpre[:] = self.Cpre


                tAll = zeros((len(tPre[1:]) + len(tPostSorted), 3))

                tAll[:, 0] = hstack((tPre[1:], tPostSorted))
                tAll[:, 1] = hstack((zeros(len(tPre[1:])), ones(len(tPostSorted))))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,len(tPostSorted))))
                tList = tAll.tolist()
                tListSorted = sorted(tList, key=lambda tList: tList[0])

                # tListSorted.append([Npres/freq,2])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD / tListSorted[-1][0]
                alphaP = tP / tListSorted[-1][0]
                # print alphaD, alphaP

                return (alphaD, alphaP)

        ###############################################################################
        # stochastic Sjoestroem 2001 protocol
        def spikePairStochasticFrequency(self,DeltaTStart,DeltaTEnd,freq,Npres):
                tStart = 0.1 # start time at 100 ms
                
                Npres = Npres*12
                timeAbove = zeros((1,2))
                
                tD = 0.
                tP = 0.
                random.seed(7)
                tPre = arange(Npres)/freq + tStart + (DeltaTStart + np.random.rand(Npres)*(DeltaTEnd-DeltaTStart))
                tPost = tPre +  (DeltaTStart + np.random.rand(Npres)*(DeltaTEnd-DeltaTStart))
                
                tAll = zeros((2*Npres,3))
                
                tAll[:,0] = hstack((tPre,tPost))
                tAll[:,1] = hstack((zeros(Npres),ones(Npres)))
                tAll[:,2] = hstack((repeat(self.Cpre, Npres), repeat(self.Cpost, Npres)))
                tList = tAll.tolist()
                
                tListSorted = sorted(tList, key=lambda tList: tList[0])
                
                tListSorted.append([Npres/freq,2])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                alphaD = tD/(float(Npres))
                alphaP = tP/(float(Npres))
                
                return (tD/float(Npres),tP/float(Npres))
            
                #return (alphaD,alphaP)
                
        ###############################################################################
        def spikePairFrequencySTPDeterministic(self, deltaT, freq, Npres, tauRec, U ):
                tStart = 0.1  # start time at 100 ms

                #Npres = Npres * 12
                #timeAbove = zeros((1, 2))

                #random.seed(7)
                tPre  = arange(Npres) / freq + tStart
                tPost = tPre + deltaT

                tAll = zeros((2 * Npres, 3))

                cpre = zeros(Npres)
                if U != 0:
                    cpre[0] = 1.
                    for i in range(1, Npres):
                        cpre[i] = 1. - (1. - (cpre[i-1] - U*cpre[i-1]))*np.exp(-1./(freq * tauRec))
                    cpre *= U*self.Cpre
                else:
                    cpre[:] = self.Cpre

                #print cpre
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npres), ones(Npres)))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,Npres)))
                tList = tAll.tolist()

                tListSorted = sorted(tList, key=lambda tList: tList[0])

                tListSorted.append([Npres / freq + tStart, 2,0])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                return (tD , tP )

        ###############################################################################
        def spikePairFrequencySTPStochastic(self, deltaT, freq, Npres, tauRec, pRelease, Nves):

                tStart = 0.1  # start time at 100 ms
                NrepetitionsStoch = 100
                q = self.Cpre / Nves

                tPre  = arange(Npres) / freq + tStart
                tPost = tPre + deltaT

                #######################################
                # stochastic STD model implementation
                # ampStoch = np.zeros((NrepetitionsStoch, len(tPre)))
                timesAbove = zeros((NrepetitionsStoch, 2))
                # tP = 0.
                for r in range(NrepetitionsStoch):
                        Vesicles = np.ones((len(tPre), Nves))
                        Release = np.random.rand(len(tPre), Nves) < pRelease
                        # VesTimes = transpose(np.tile(tPre, (Nves, 1)))
                        for i in range(len(tPre)):
                                nRel = sum(Release[i])
                                releaseSites = argwhere(Release[i] == True)
                                emptyTimes = random.exponential(tauRec, nRel)
                                for n in range(nRel):
                                        if not Vesicles[i, releaseSites[n]] == 0.:
                                                recoveryMask = ((tPre - tPre[i]) < emptyTimes[n]) & (
                                                                (tPre - tPre[i]) > 0)
                                                Vesicles[recoveryMask, releaseSites[n]] = 0
                        # print 'after amp det'
                        cpreStoch = q * (sum(Vesicles * Release, axis=1))

                        ######################################
                        tAll = zeros((len(tPre) + len(tPost), 3))

                        tAll[:, 0] = hstack((tPre, tPost))
                        tAll[:, 1] = hstack((zeros(len(tPre)), ones(len(tPost))))
                        tAll[:, 2] = hstack((cpreStoch, repeat(self.Cpost, len(tPost))))
                        tList = tAll.tolist()
                        tListSorted = sorted(tList, key=lambda tList: tList[0])

                        tListSorted.append([Npres / freq + tStart, 2, 0])

                        # tListSorted.append([Npres/freq,2])

                        ###########################################################
                        # event-based integration
                        (timesAbove[r, 0], timesAbove[r, 1]) = self.eventBasedIntegration(tListSorted)
                        # print r
                return (average(timesAbove[:,0]) , average(timesAbove[:,1]) )


        ###############################################################################
        # (timeD,timeP) = tat.spikePairFrequencyNonlinear(DeltaTStart,DeltaTEnd,D,frequency)
        def spikePairFrequencyTMM(self, deltaT, freq, Npres, U, tauFac, tauDep ):
                tStart = 0.1  # start time at 100 ms

                #Npres = Npres * 12
                #timeAbove = zeros((1, 2))

                tD = 0.
                tP = 0.
                #random.seed(7)
                tPre  = arange(Npres) / freq + tStart
                tPost = tPre + deltaT

                tAll = zeros((2 * Npres, 3))

                cpre = zeros(Npres)
                u   = zeros(Npres)
                x   = zeros(Npres)
                if U != 0:
                        u[0] = U
                        x[0] = 1.
                        cpre[0] = u[0] * x[0]
                        for i in range(1, len(tPre[1:])):
                                u[i] = u[i - 1] * np.exp(-1./(freq*tauFac)) + U * (1. - u[i - 1] * np.exp(-1./(freq*tauFac)))
                                x[i] = 1. - (1. - (x[i - 1] - u[i - 1] * x[i - 1])) * np.exp(
                                        -1./(freq*tauDep))
                                cpre[i] = u[i] * x[i]
                        cpre *= self.Cpre
                        #cpre[0] = 1.
                        #for i in range(1, Npres):
                        #    cpre[i] = 1. - (1. - (cpre[i-1] - U*cpre[i-1]))*np.exp(-1./(freq * tauRec))
                        #cpre *= U*self.Cpre
                else:
                        cpre[:] = self.Cpre

                #print cpre
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npres), ones(Npres)))
                tAll[:, 2] = hstack((cpre,repeat(self.Cpost,Npres)))
                tList = tAll.tolist()

                tListSorted = sorted(tList, key=lambda tList: tList[0])

                tListSorted.append([Npres / freq + tStart, 2,0])

                ###########################################################
                # event-based integration
                (tD, tP) = self.eventBasedIntegration(tListSorted)

                return (tD , tP )

        ###############################################################################
        # (timeD,timeP) = tat.spikePairFrequencyNonlinear(DeltaTStart,DeltaTEnd,D,frequency)
        def spikePairStochasticFrequencySTPDeterministic(self, DeltaTStart, DeltaTEnd, freq, Npres, tauRec, U):
            tStart = 0.1  # start time at 100 ms

            # Npres = Npres * 12
            # timeAbove = zeros((1, 2))
            Naverages = 100.
            tD = 0.
            tP = 0.
            for i in arange(Naverages):
                # random.seed(7)
                tPre = arange(Npres) / freq + tStart + (DeltaTStart + rand(Npres) * (DeltaTEnd - DeltaTStart))
                tPost = tPre + (DeltaTStart + rand(Npres) * (DeltaTEnd - DeltaTStart))

                # tPre  = arange(Npres) / freq + tStart
                # tPost = tPre + deltaT

                tAll = zeros((2 * Npres, 3))

                cpre = zeros(Npres)
                cpre[0] = self.Cpre * U
                for i in range(1, Npres):
                    cpre[i] = cpre[i - 1] * (1. - U * np.exp(-(tPre[i] - tPre[i - 1]) / (tauRec)))
                tAll[:, 0] = hstack((tPre, tPost))
                tAll[:, 1] = hstack((zeros(Npres), ones(Npres)))
                tAll[:, 2] = hstack((cpre, repeat(self.Cpost, Npres)))
                tList = tAll.tolist()

                tListSorted = sorted(tList, key=lambda tList: tList[0])

                tListSorted.append([Npres / freq + tStart, 2])

                ###########################################################
                # event-based integration
                (tDTemp, tPTemp) = self.eventBasedIntegration(tListSorted)

                tD += tDTemp
                tP += tPTemp

            return (tD / Naverages, tP / Naverages)



        ###############################################################################
        # (timeD,timeP) = tat.spikePairFrequencyNonlinear(DeltaTStart,DeltaTEnd,D,frequency)
        def eventBasedIntegration(self, tListSorted):

                # event-based integration
                ca = []
                # CaTotal, CaPre, CaPost, time
                ca.append([0.,0.,0.,0.])
                pre = 0
                post = 0
                tD = 0.
                tP = 0.
                for i in tListSorted:
                        #
                        caTotOld    = ca[-1][0]
                        caPreOld    = ca[-1][1]
                        caPostOld   = ca[-1][2]
                        tOld        = ca[-1][3]
                        #caTotTemp  = caTotOld*np.exp(-(i[0]-tOld)/self.tauCa)
                        caPreTemp  = caPreOld*np.exp(-(i[0]-tOld)/self.tauCa)
                        caPostTemp = caPostOld*np.exp(-(i[0]-tOld)/self.tauCa)
                        caTotTemp  = caPreTemp + caPostTemp
                        if caTotOld > self.thetaD:
                                if caTotTemp > self.thetaD:
                                        tD += i[0]-tOld
                                else:
                                        tD += (self.tauCa)*np.log(caTotOld/self.thetaD)
                        if caTotOld > self.thetaP:
                                if caTotTemp > self.thetaP:
                                        tP += i[0]-tOld
                                else:
                                        tP += (self.tauCa)*np.log(caTotOld/self.thetaP)
                        # postsynaptic spike
                        if i[1] == 1:
                                caPostTemp += i[2] + self.eta*caPreTemp
                                post+=1
                        # presynaptic spike
                        if i[1] == 0:
                                caPreTemp += i[2] #self.Cpre
                                pre+=1
                        caTotTemp = caPreTemp + caPostTemp
                        ca.append([caTotTemp,caPreTemp,caPostTemp,i[0]])
                        #
                        #pdb.set_trace()

                return (tD, tP)
