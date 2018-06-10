#!/usr/bin/env python
from mapHrEstimator import *


#
# Global function
#

class Tracker:

    def __init__(self, start, alpha=.01, beta=0, deltaFreqState = np.float(0), time=-1000,
                 maxChange = .5, boundHi=205, boundLo=40, maxDeltaT=3000):

        self.freqState = np.float(start)
        self.deltaFreqState = deltaFreqState
        self.time = time

        self.boundHi = boundHi
        self.boundLo = boundLo

        self.peakHist = []
        self.freq = []
        self.deltaFreq = []
        self.timeHist = []
        self.drHist = []

        self.alpha = alpha
        self.beta = beta
        self.maxChange = maxChange
        self.maxDeltaT = maxDeltaT

    def update(self, time, peak, dynamicRange=None, maxRes=20):

        deltaT = (time - self.time)

        if deltaT > self.maxDeltaT:
            deltaT = self.maxDeltaT

        #Convert into seconds
        deltaT = deltaT/1000


        #todo - why do we need this???
        if deltaT <= 0.0:
            print("Negative DeltaT")
            return 0


        self.time = time
        self.timeHist.append(self.time)
        self.drHist.append(dynamicRange)

        if peak == -1:
            self.setInvalidHR(invalidHRHold=invalidHRHold)
            return 0

        if peak is None:
            print("No Peak Passed to tracker")
            self.peakHist.append(0)
        else:
            self.peakHist.append(peak)

        if peak is not None:
            if peak < self.boundLo or peak > self.boundHi:
                peak = self.freqState
                self.deltaFreqState = 0
        else:
            self.deltaFreqState = 0



        if self.deltaFreqState > .5:
            self.deltaFreqState = .5

        if self.deltaFreqState < -.5:
            self.deltaFreqState = -.5


        # Kludge: Setting deltaFreqState to zero thus eliminated the beta part of the filter
        self.deltaFreqState = 0

        self.freqState += deltaT*self.deltaFreqState

        if peak is not None:

            residual = peak - self.freqState


            alpha = self.alpha
            beta = self.beta

            if np.abs(residual) > maxRes:
                residual = np.sign(residual)*maxRes

            #update the state
            self.freqState += alpha*residual
            self.deltaFreqState += (beta/deltaT)*residual

        if self.freqState < self.boundLo:
            self.freqState = self.boundLo
            self.deltaFreqState = 0

        elif self.freqState > self.boundHi:
            self.freqState = self.boundHi
            self.deltaFreqState = 0


        self.freq.append(self.freqState)
        self.deltaFreq.append(self.deltaFreqState)


        return 0


    def setInvalidHR(self, invalidHRHold=False):
        self.deltaFreqState = 0
        self.peakHist.append(0)
        if invalidHRHold:
            self.freq.append(self.freqState)  # hold prevHR during HR is invalid
        else:
            self.freq.append(-1)  # do not hold prevHR, output -1

        self.deltaFreq.append(self.deltaFreqState)

        return 0



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('fname', help='data log file', default=None)
    parser.add_argument('--truthFile', help='heart strap data', default=None)
    parser.add_argument('--alphaFile', help='heart strap data from mio alpha', default=None)
    parser.add_argument('--noPlots',help='show the plots or not',default=False,action='store_true')
    parser.add_argument('--out', help='output filename', default='foo.csv')
    args = parser.parse_args()


    optTracks = []
    accTracks = []

    #header_def = [ ('time',float),('opt0',float),('opt1',float),('opt2',float),('acc0',float),('acc1',float),('acc2',float) ]
    d = np.genfromtxt(args.fname,delimiter=',')
    peaks = d[:,1:4]
    accPeaks = d[:,4:7]
    time = d[:,0]

    startVals = [70]
    for ind in np.arange(0,len(startVals)):
        optTracks.append(Tracker(startVals[ind],maxChange=5))

    startVals = [100]
    for ind in np.arange(0,len(startVals)):
        accTracks.append(Tracker(startVals[ind],alpha=.1,beta=.25))


    for ind in np.arange(0,peaks.shape[0]):



        for peakInd in np.arange(0,peaks.shape[1]):

            for accInd in np.arange(0,accPeaks.shape[1]):
                if (np.abs(peaks[ind,0] - peaks[ind,1]) < 20):
                    if np.abs(peaks[ind,peakInd]-accPeaks[ind,accInd]) < 5.0:
                        peaks[ind,peakInd] = np.min(peaks[ind,:])


        #update the accelerometer tracks
        #for each track find the closest peak
        for track in np.arange(0,len(accTracks)):
            accTracks[track].update(time[ind],accPeaks[ind,track])
        '''
        #for each track find the closest peak
        for track in accTracks:

            res = np.zeros((accPeaks.shape[1],))
            for peakInd in np.arange(0,accPeaks.shape[1]):
                res[peakInd] =  np.abs(accPeaks[ind,peakInd] - track.freqState)

            closest = np.argmin(res)
            track.update(time[ind],accPeaks[ind,closest])
        '''


        #for each track find the closest peak
        for track in optTracks:

            res = np.zeros((peaks.shape[1],))
            weight=np.array([1.0,1.0,1.0])
            for peakInd in np.arange(0,peaks.shape[1]):
                if peaks[ind,peakInd] > 90:
                    res[peakInd] =  weight[peakInd]*np.abs(peaks[ind,peakInd] - track.freqState)

            closest = np.argmin(res)
            track.update(time[ind],peaks[ind,closest])


    pl.figure()
    for ind in np.arange(0,peaks.shape[1]):
         pl.plot(time[:],peaks[:,ind],'+')
    pl.grid(True)

    #pl.figure()
    #todo - interpolate truth heart rate onto measured heart rate
    if args.truthFile is not None:
        hrTruth = np.genfromtxt(args.truthFile,skiprows=3,delimiter=',');
        tTrue=hrTruth[:,0]-hrTruth[1,0]
        tTrue /= 1000
        pl.plot(tTrue,hrTruth[:,1],'g')


    for track in optTracks:
        pl.plot(track.timeHist,track.freq,'--')
    pl.grid(True)


    pl.figure()
    for ind in np.arange(0,accPeaks.shape[1]):
         pl.plot(time[:],accPeaks[:,ind],'+')

    for track in accTracks:
        pl.plot(track.timeHist,track.freq,'--')
    pl.grid(True)


    pl.figure()
    if args.truthFile is not None:
        pl.plot(tTrue,hrTruth[:,1],'g')
    for track in optTracks:
        pl.plot(track.timeHist,track.freq,'--')
    pl.grid(True)

    pl.figure()
    pl.plot(optTracks[0].residual)


    pl.show()

