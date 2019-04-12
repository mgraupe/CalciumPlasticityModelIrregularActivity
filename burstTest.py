import numpy as np
import pdb

def separateBursts(spikeTimes, burstInterval):

    # extract bursts ###############################
    spikeTimesA = np.asarray(spikeTimes)
    diffSpikeTimes = np.diff(spikeTimesA)
    boolBursts = diffSpikeTimes < burstInterval
    # displace bool array to include also last spikes in bursts
    boolBurstsAll = np.concatenate((boolBursts,np.array([False]))) | np.concatenate((np.array([False]),boolBursts))
    onlyBursts = spikeTimesA[boolBurstsAll]
    # extract single spikes ########################
    boolBurstsAllInt = np.array(boolBursts,dtype=int)
    diffBurst = np.diff(boolBurstsAllInt)
    noBurstsBool = (np.invert(boolBurstsAll)) | (np.concatenate((np.array([1]),diffBurst,np.array([0])))==1)
    noBursts = spikeTimesA[noBurstsBool]

    # count bursts and number of spikes in bursts
    diffB = onlyBursts[1:] - onlyBursts[:-1]
    numberOfBursts = np.sum(diffB > burstInterval) + 1  # number of bursts is given by intervals larger than the 'burstInterval'
    boolBursts = diffB < burstInterval
    boolBursts = np.concatenate((np.array([False]), boolBursts))
    boolBursts = np.concatenate((np.diff(boolBursts), np.array([True])))
    startStopBursts = np.arange(len(boolBursts))[boolBursts]
    # [i for i,(m,n) in enumerate(zip([2]+boolBursts,boolBursts+[2])) if m!=n]
    # startStopBursts = np.asarray(startStopBursts)
    numberSpikesInBursts = (startStopBursts[1::2] - startStopBursts[:-1:2]) + 1
    # pdb.set_trace()
    return (onlyBursts, noBursts, numberOfBursts, numberSpikesInBursts)


nSpikes = 1000
preRate = 1.

tPre = []
tPre.append(0)

for i in range(int(nSpikes)):
    tPre.append(tPre[-1] + np.random.exponential(1. / preRate))

(tPreBursts,_,numberOfBursts,numberSpikesInBursts)  = separateBursts(tPre[1:],0.15)

