#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as pl
import signal as sig
import scipy.signal as signal
import sensorLog as log
from collections import OrderedDict
import iirFiltDesign as iir

def signal_handler(signal, frame):
    sys.exit(0)
sig.signal(sig.SIGINT, sig.SIG_DFL)

class PowerSpec:

    def __init__(self, fftLen=8192,sampleRate=100,winLen=1000,peakHi_fft=205,peakLo_fft=17,peakHi_hps=800,peakLo_hps=17,
                 numHarmonics=4,numPeaks=None, plotBPMHigh=300, plotBPMLow=10, useWelch=False, use_threshold=True, dBcutoff=20,
                 dBcutoff_hps=10, label=""):

        self.label= label
        self.fftLen = fftLen
        self.sampleRate = sampleRate
        self.winLen = winLen
        self.freq = np.fft.rfftfreq(fftLen, d=1/sampleRate)
        self.window = signal.blackman(self.winLen, sym=True)

        self.numHarmonics = numHarmonics
        self.numPeaks = numPeaks

        # Peak boundaries for fft peak search
        self.peakHi_fft = peakHi_fft
        self.peakLo_fft = peakLo_fft
        self.bpmRange_fft, self.bpmBins_fft = self.freq2bpmRange(self.peakLo_fft, self.peakHi_fft)

        # Peak boundaries for hps peak search
        self.peakHi_hps = peakHi_hps
        self.peakLo_hps = peakLo_hps
        self.bpmRange_hps, self.bpmBins_hps = self.freq2bpmRange(self.peakLo_fft, self.peakHi_fft)

        #Plotting boundaries
        self.plotBPMHigh = plotBPMHigh
        self.plotBPMLow = plotBPMLow
        self.plot_bpmRange, self.plot_bpmBins = self.freq2bpmRange(self.plotBPMLow, self.plotBPMHigh)

        self.use_threshold = use_threshold
        self.dBcutoff = dBcutoff
        self.dBcutoff_hps = dBcutoff_hps

        self.useWelch = useWelch

        self.fft_har_hist = []
        self.hps_har_hist = []
        self.har_hist = []

        self.harmonic = False
        self.fft_harmonic = False
        self.hps_harmonic = False

        self.hps_fundamental = 0
        self.fft_fundamental = 0

        self.tolerance = 5

        self.topPeak = 0
        self.notchedFreq = []


    def freq2bpmRange(self, freqLo, freqHi):
        '''
        Given HR values ranging from freqLo to freqHi, this function returns the bins for this range.

        :param freqLo: Lower frequency bound
        :param freqHi: Upper frequency bound
        :return:
            bpmBins: index of bpm bins
            bpmRange: frequency value of each bpm bin
        '''

        hiBin = self.bpm2bin(freqHi)
        loBin = self.bpm2bin(freqLo)
        bpmBins = np.arange(loBin,hiBin)
        bpmRange = (bpmBins * self.sampleRate * 60.0 / self.fftLen)

        return bpmRange, bpmBins


    def bin2bpm(self,bin):
        '''
        This function converts a bin index to bpm
        :param bin: input bin index
        :return: HR in bpm
        '''

        return int(bin*self.sampleRate*60.0/self.fftLen)


    def bpm2bin(self,bpm):
        '''
        This function converts a bpm value to a bin index number
        :param bin: in bpm
        :return: HR input bin index
        '''

        return int(bpm*self.fftLen/(self.sampleRate*60)+0.5)


    def getPowerSpec(self,x):
        '''
        This function generates the FFT of a given signal x, and then finds all the relevant peaks
        It populates specPeaks and specVals
        :param x: Time domain signal
        :return:
            self.X: input signal FFT
            self.specVals: magnitude of relevant peaks found (ordered in descending order)
            self.specPeaks: relevant frequency peaks found  (match magnitude)
        '''

        self.x = np.copy(x)

        if self.useWelch:
            f, self.X = signal.welch(self.x, self.sampleRate, nperseg=1024,window='blackman')
        else:
            self.X = 2*np.abs(np.fft.rfft(x*self.window, n=self.fftLen))
            # Multiplying by 2 is a scaling factor because we are only looking at the positive side of the spectrum

        self.specPeaks, self.specVals = self.getLocalPeaks(self.X, self.dBcutoff,peakLo=self.peakLo_fft,
                                                           peakHi=self.peakHi_fft, searchBpmBins=self.bpmBins_fft,
                                                           numPeaks=self.numPeaks)

        return



    def getLocalPeaks(self,x,dBcutoff,peakLo,peakHi,searchBpmBins,numPeaks=None):

        '''
        This function finds the top peaks of the input x, which can be the HPS or FFT of a signal.
        The peaks are searched in the area of interest (searchBpmBins)
        The peaks are found within the given threshold (dBcufoff) and boundaries (peakLo and peakHi).
        The return is ordered by decreasing order of magnitude.

        :param x: HPS or FFT of a signal
        :param dBcutoff: dB threshold below which peaks are ignored
        :param peakLo: lower frequency bound
        :param peakHi: higher frequency bound
        :param searchBpmBins: area of search
        :param numPeaks: optional argument, maximum number of peaks to be returned
                        if nothing is passed, it returns all peaks found
        :return:
                peakFreqs: bpm frequency of FFT peaks found give the constraints
                peakMag: magintude of FFT peaks found
        '''
        
        # MAX Values
        x_bpmRange = x[searchBpmBins]
        maxPeak = x_bpmRange.argmax()
        maxPeakVal = x_bpmRange[maxPeak]
        maxPeak_dB= 20*np.log10(maxPeakVal)


        # find the relative peaks in the spectrum
        peakLocations = signal.argrelmax(x)
        peakVals = x[peakLocations]

        # peakMap stores the frequency to peak magnitude correspondence
        peakMap = dict(zip(self.freq[peakLocations], peakVals))

        # There is a chance that the maxPeak itself is un-physiological
        # In that case, it will get eliminated below.
        badCandidates = set()
        for frequency, peakVal in peakMap.items():
            # IF the value does not meet our condition, note it
            if not self.isCandidate(frequency, peakVal, maxPeak_dB, dBcutoff, peakLo, peakHi):
                badCandidates.add(frequency)

        # Remove the bad candidates
        for frequency in badCandidates:
            peakMap.pop(frequency)

        #
        # Sort by dictionary value, frequencies sorted by their strengths
        #
        ordered = OrderedDict(sorted(peakMap.items(), key=lambda d: d[1], reverse=True))
        counter = 0

        if numPeaks is None:
            numPeaks = len(ordered)

        peakFreqs= np.zeros((numPeaks,))
        peakMags = np.zeros((numPeaks,))

        for freq in ordered:
            # We want to pick the top three frequencies in each channel
            if counter >= numPeaks:
                break
            # orderedFreqsByIntensity[freq]) stores the intensity of the peak in dB
            peakFreqs[counter] = 60*freq
            peakMags[counter] = ordered[freq]
            counter += 1

        return peakFreqs,peakMags


    def isCandidate(self, freq, peakValue, maxPeakDB, dBcutoff, peakLo, peakHi):
        '''
        This function tests if a given freq value falls within the boundaries and is within
        dB threshold
        :param freq: input frequency to be tested
        :param peakValue: amplitude of input frequency
        :param maxPeakDB: maximum peak value of search area
        :param dBcutoff: dB threshold below which freq is ignored
        :param peakLo: lower frequency bound
        :param peakHi: upper frequency bound
        :return:
            Returns True if all conditions are met (magnitude is within dB threshold,
            and freq is within boundaries). False otherwise
        '''

        if peakLo <= (freq*60) <= peakHi:

            if self.use_threshold:

                peakDB = 20*np.log10(peakValue)
                diff = np.abs(maxPeakDB - peakDB)
                if  diff <= dBcutoff:
                    return True
                else:
                    return False
            else:
                return True
        else:

            return False


    def getHPS_C(self, input):
        '''
        This is a python implementation of the C function for getHPS.
        It is here for debug and comparison purposes.
        :param input: FFT of a signal
        :return:
            out: HPS of a signal
        '''

        out = np.zeros((len(input)))

        for i in range(0,int(self.fftLen/(2*self.numHarmonics))):
            tmp = 1
            for j in range(1,self.numHarmonics+1):

                index = i*j

                if (index>self.fftLen/2):
                    break;

                tmp *= input[index]

            out[i] = tmp

        return out


    def getHPS(self):
        '''
        This function generates the HPS of a signal self.X
        :return:
            self.HPS : HPS based on self.X (fft of input signal self.x)
        '''
        self.HPS = np.array(self.X)

        for i in np.arange(1,self.numHarmonics):
            temp = np.zeros((self.X.shape[0],))
            length = len(self.X[0:len(self.X):i+1])
            temp[0:length] = self.X[0:len(self.X):i+1]
            self.HPS *= temp

        return


    def getHarmProdSpec(self):

        '''
        This function generates the Harmonic Product Spectrum of a signal (self.X), and then finds
        the relevant peaks and its amplitudes. Lastly, it determines if the FFT and HPS of the signal
		are harmonic; and if so, it also generates a fundamental frequency.

		:return:
            self.HPS : HPS based on self.X (fft of input signal self.x)
            self.hpsFreqPeaks : bpm frequency of HPS peaks found give the constraints
            self.hpsPeakAmplitudes : magnitude of HPS peaks found
        '''

        self.getHPS()

        self.hpsFreqPeaks, self.hpsPeakAmplitudes = self.getLocalPeaks(self.HPS,dBcutoff=self.dBcutoff_hps,peakLo=self.peakLo_hps,
                                                         peakHi=self.peakHi_hps, searchBpmBins=self.bpmBins_hps,numPeaks=self.numPeaks)
        self.isHarmonic()

        return



    def isHarmonic(self):
        '''
        This function determines whether or not a signal is harmonic
        This is done by first finding if a signal seems to be harmonic from its FFT peaks
        Then, we check if it seems harmonic from its HPS peaks.
        Both checks generate a suspected fundamental frequency in bpm.
        Finally, if the signal seems harmonic in both cases, we check if the suspected fundamentals
        match or are very close to each other or are multiple of each other.

        :return: True if all checks are passed; therefore signal is harmonic, or False otherwise.
            self.fft_fundamental = 0
            self.fft_harmonic = False
            self.fft_fundamental : Fundamental frequency in FFT
            self.fft_harmonic : True is FFT seems harmonic, false otherwise
            self.hps_fundamental : Fundamental frequency in HPS
            self.hps_harmonic : True is HPS seems harmonic, false otherwise
            self.harmonic : True is both fft and hps harmonic, and both suspected fundamentals are near
                            each other or are multiples of each other within tolerance
            self.har_hist : List of previous True/False decision on self.harmonic
            self.fft_har_hist : List of previous True/False decision on self.fft_harmonic
            self.hps_har_hist : List of previous True/False decision on self.hps_harmonic
        '''

        self.is_FFT_Harmonic()
        self.is_HPS_Harmonic()

        self.harmonic = False
        self.fundamental = 0

        if self.fft_harmonic and self.hps_harmonic:
            if np.abs(self.hps_fundamental-self.fft_fundamental)<= self.tolerance:
                self.fundamental = self.fft_fundamental
                self.harmonic=True

        self.har_hist.append(self.harmonic)
        self.hps_har_hist.append(self.hps_harmonic)
        self.fft_har_hist.append(self.fft_harmonic)

        return self.harmonic



    def is_HPS_Harmonic(self):
        '''
        is_HPS_Harmonic decides if the HPS of a signal is harmonic and then finds a possible fundamental for such signal

        It does so by checking if there is at least one hps peak but no more than two.
        If zero peaks are found, the function returns immediately.
        If there are more thank two peaks, the function also returns immediately.

        Then, we proceed to check that the chosen peak is not just a small blip in a decaying exponential
        This is done by checking if the top peak is the greatest value in the area of interested.
        If not, then we check if the top peak is very close (<=5bpm) to the lower frequency bound.
        If this check is False, the function also returns.

        If there is exactly one peak, this peak is chosen as the HPS peak.

        If there are exactly two peaks, and additional check is performed:
        In general (99% of time), the tallest peak is the fundamental; however, in some occasions the 
		fundamental isn't the tallest peak, but the second tallest. That's why it is necessary to check if the 
		second tallest peak is the fundamental, and the tallest is a harmonic (a multiple of the lower peak).

        If 0 peaks: hps_fundamental=0; hps_harmonic=False
        If >2 peaks: hps_fundamental=0; hps_harmonic=False
        Check if peak is valid (not a decaying exponential)
        If 1 peak: Peak is fundamental and hps_harmonic=True
        If 2 peaks: Check if tallest, check if the lowest frequency peak is the tallest peak
                    If not, check that if the lower peak of lower frequency is the fundamental and
                    tallest peak of higher frequency is the harmonic.

        :returns:
            self.fft_fundamental : Suspected HPS fundamental
            self.fft_harmonic : True is HPS seems harmonic, false otherwise
        '''

        self.hps_harmonic = False

        #Remove Zeros from hpsPeaks
        hps_peaks = [i for i in self.hpsFreqPeaks if i != 0]

        self.hps_fundamental = 0

        # If there are no peaks, return
        if len(hps_peaks)==0:
             return


        # Check that there is at least one HPS peak, but no more than two.
        # This ensures that HPS isn't too cluttered : i.e. there a lots of peak in a region of interest.
        # HPS should ideally have only a single strong peak
        if len(hps_peaks)>2:
            return


        # Sometime the HPS will look like a decaying exponential, and we may pick up a zig in lower frequency bound
        # in this case, the HPS peak is not interesting. We exit out of the function if the peak found is not the highest
        # value of our area of interest and if it is too close to the lower frenquency bound.
        if (self.hpsPeakAmplitudes[0] < np.max(self.HPS[self.bpmBins_hps])) and (self.hpsFreqPeaks[0] <= self.peakLo_hps+self.tolerance):
            return

        # At this point we have one or two peak, so we assume HPS is harmonic.
        self.hps_harmonic = True

        # Return if there is a single peak
        if len(hps_peaks) == 1:
            self.hps_fundamental = self.hpsFreqPeaks[0]
            return


        # At this point, we have exactly two hps peaks. Test if the lower frequency peak has a lower amplitude
        if self.hpsFreqPeaks[0] > self.hpsFreqPeaks[1]:
            # Note: the spectral peaks have been sorted by peak amplitude (descending order)
            
            # If the higher frequency peak is a multiple of the lower frequency peak, select the lower frequency peak.
            residual = self.hpsFreqPeaks[0] % self.hpsFreqPeaks[1]
            if residual <= self.tolerance or residual >= (self.hpsFreqPeaks[1] - self.tolerance):
                self.hps_fundamental = self.hpsFreqPeaks[1]
            # If not a multiple, select the higher frequency peak.
            else:
                self.hps_fundamental = self.hpsFreqPeaks[0]
        else:
            self.hps_fundamental = self.hpsFreqPeaks[0]

        return


    def is_FFT_Harmonic(self):

        '''
        is_FFT_Harmonic determines if a signal shows a harmonic FFT
        is_FFT_Harmonic finds a possible fundamental frequency from specPeaks, and checks if all other peaks
        in specPeaks are harmonics (or multiples) of it.
        The fundamental is smallest peak value in specPeaks or the smallest delta between adjacent peaks, whichever
        one is smaller.
        Once the fundamental is found, we check if the rest of the peaks are multiples of it, and are within tolerance.
        If any single other peak is not a multiple of it, we assume the signal is not harmonic
        :return
            self.fft_fundamental : Suspected FFT fundamental
            self.fft_harmonic : True is FFT seems harmonic, false otherwise
        '''

        self.fft_fundamental = 0
        self.fft_harmonic = False

        #Remove Zeros from fft_peaks and sort
        freq_peaks = [i for i in self.specPeaks if i !=0]
        sorted_freq_peaks = sorted(freq_peaks)

        # If there is only one freq_peak, assume it is harmonic
        if len(sorted_freq_peaks) == 1:
            # Note: The fft_fundamental_value is bounded by self.peakLo_fft and self.peakHi_fft
            self.fft_fundamental = sorted_freq_peaks[0]
            self.fft_harmonic = True

        elif len(sorted_freq_peaks) > 1:

            #Get possible fundamental: Smallest freq peak in the list
            temp_fft_fundamental = np.min(sorted_freq_peaks)

            if len(sorted_freq_peaks) == 2:
                #Fundamental could also be the interval between two peaks:
                temp_fft_fundamental = np.min([np.diff(sorted_freq_peaks), temp_fft_fundamental])

            else:
                possible_fundamentals_deltas = np.diff(sorted_freq_peaks)
                #The smallest delta in frequency peaks must be the fundamental:
                temp_fft_fundamental = np.min([np.min(possible_fundamentals_deltas), temp_fft_fundamental])

            # Check that fundamental is larger that lowest acceptable threshold:
            if temp_fft_fundamental >= self.peakLo_fft:

                #Peak may not line up perfectly, so check that ALL residuals are small enough to consider them harmonics
                residuals = np.abs([i%temp_fft_fundamental for i in sorted_freq_peaks])

                #harmonic_test gets populated with 0 if peak is harmonic or 1 is it isn't
                harmonic_test = np.ones((len(residuals)))

                for i in np.arange(len(residuals)):
                    if residuals[i] <= self.tolerance or (temp_fft_fundamental-self.tolerance) <= residuals[i]:
                        harmonic_test[i] = 0

                # Check that all peaks were harmonics, if so return the fun
                if np.sum(harmonic_test) == 0:
                    self.fft_harmonic = True
                    self.fft_fundamental = temp_fft_fundamental

        return


    def getTopPeak(self, comparisonPeaks, dBDistance, bpmDistanceToAccPeaks):

        # This function checks if the top optical peak is much larger than the next highest peak.
        # Also, the top peak must not exist in the accelerometer
        # Returns 0 if conditions were not met, or the top-peak otherwise.

        self.topPeak = 0

        if (20*np.log10(self.specVals[0])-20*np.log10(self.specVals[1])) >= dBDistance:

            bpmDistance = [np.abs(i-self.specPeaks[0]) for i in comparisonPeaks]

            if np.min(bpmDistance) >= bpmDistanceToAccPeaks:

                self.topPeak = self.specPeaks[0]

        return self.topPeak


    def plotdB(self, title=""):

        title = title + self.label + ' in dB'
        f1, ax1 = pl.subplots(2, sharex=True)
        f1.suptitle(title)

        ax1[0].plot(self.plot_bpmRange, 20*np.log10(self.X[self.plot_bpmBins]), 'b')
        ax1[0].set_ylabel('Power Spectrum dB')
        ax1[1].set_xlabel('bpm')
        ax1[0].grid()

        ax1[1].plot(self.plot_bpmRange, 20*np.log10(self.HPS[self.plot_bpmBins]), 'b')
        ax1[1].set_ylabel('Harmonic Product Spectrum dB')
        ax1[1].set_xlabel('bpm')
        ax1[1].grid()

        return f1


    def plot(self, title="", ax=None):

        title = title + self.label

        f1, ax = pl.subplots(2, sharex=True)
        f1.suptitle(title)

        self.plotSpectrum(ax=ax[0])
        self.plotHarmonicSpectrum(ax=ax[1])

        return f1


    def plotTD(self, title="", ax=None, sharex=None, color='b'):
        title = title + self.label + ' TD Signal'
        if ax is None :
            fig, ax = pl.subplots()
        time_sec = np.arange(0,1000)/100
        if sharex is None:
            ax.plot(time_sec, self.x, color=color)
        else:
            ax.plot(time_sec, self.x, sharex=sharex, color=color)
        ln = ax.get_lines()[-1]
        ax.set_title(title)
        ax.set_xlabel('Time (sec)')
        ax.grid()

        return ln


    def plotHarmonicSpectrum(self, ax, title="", sharex=None, color='b'):
        if sharex is None:
            ax.plot(self.plot_bpmRange, 20*np.log10(self.HPS[self.plot_bpmBins]), color=color)
        else:
            ax.plot(self.plot_bpmRange, 20*np.log10(self.HPS[self.plot_bpmBins]), color=color, sharex=sharex)

        ax.plot(self.hpsFreqPeaks, 20*np.log10(self.hpsPeakAmplitudes), '.', ms=8, color=color)
        ax.axvline(self.hps_fundamental, color=color)
        ax.legend([self.label+' HPS (dB)'], fontsize=8, loc='best')
        ax.set_xlabel('bpm')
        return ax


    def plotSpectrum(self, ax, title="", sharex=None, truth=None, color='b'):
        if sharex is None:
            ax.plot(self.plot_bpmRange, 20*np.log10(self.X[self.plot_bpmBins]), color=color)
        else:
            ax.plot(self.plot_bpmRange, 20*np.log10(self.X[self.plot_bpmBins]), color=color, sharex=sharex)
        ln = ax.get_lines()[-1]

        ax.plot(self.specPeaks, 20*np.log10(self.specVals), '.', ms=8, color=color)
        ax.axvline(self.fft_fundamental, color=color)
        if truth is not None:
            ax.axvline(truth, color='r')
        ax.legend([self.label+' FFT (dB)'], fontsize=8, loc='best')
        return ln

    def plotNotchedFreq(self, ax, color='b'):
        '''
        plot the notched frequency onto specified axes
        :param ax:
        :param color:
        :return:
        '''
        if self.fft_harmonic and self.hps_harmonic:
            for f in self.notchedFreq:
                ax.axvline(f, color=color, ls=':', lw=1)

    def plotHarmonicHist(self, title=""):

        title = title + self.label + ' Harmonic History'
        f1, ax1 = pl.subplots(3,sharex=True)

        f1.suptitle(title)

        ax1[0].plot(self.fft_har_hist)
        ax1[0].set_ylim([-0.05,1.05])
        ax1[0].grid()
        ax1[0].set_xlabel('Time')
        ax1[0].legend(['FFT HAR'])

        ax1[1].plot(self.hps_har_hist)
        ax1[1].set_ylim([-0.05,1.05])
        ax1[1].grid()
        ax1[1].set_xlabel('Time')
        ax1[1].legend(['HPS HAR'])

        ax1[2].plot(self.har_hist)
        ax1[2].set_ylim([-0.05,1.05])
        ax1[2].grid()
        ax1[2].set_xlabel('Time')
        ax1[2].legend(['HAR'])

        return f1

    def getNotchedFreq(self, maxNotch):
        '''
        compute the notched frequencies from hps_fundamental
        :param maxNotch:
        :return:
        '''
        self.notchedFreq = []
        tmp = self.hps_fundamental
        while tmp < maxNotch and len(self.notchedFreq)<10:
            self.notchedFreq.append(tmp)
            tmp += self.hps_fundamental
        return True

    def getSpecStrength2(self, refHR, tolHR=6, fmt=None):
        '''
        compute the spectral amplitude of the signal around the reference point (5bpm diff). If there is no local peak
        in the area of interest, return invalid.
        :param ref: array [nx1] reference point (bmp).
        :param fmt: format for screen print out
        :return:
        + pk: nPk x 4 array; each row contains [ref point, minDis, pk, pkVal]
            o ref point(bpm:
            o minDis(bpm): distance from the spectral peak to ref point
            o pk(bpm): location of the spectral peak
            o pkVal(dB): spectral amplitude of the spectral peak
        + strOut
        '''
        pk = []
        strOut = ''
        i = -1
        for iref in refHR:
            i += 1
            bpmRange, bpmBin = self.freq2bpmRange(iref-tolHR, iref+tolHR)
            Px = self.X[bpmBin]
            mxIdx = signal.argrelmax(Px)[0] # location of all local peaks
            if len(mxIdx)>0: # if there are local peaks
                mnDis = tolHR
                mnDisIdx = 100
                for idxTmp in mxIdx:  # get the local peak that closest to trueHR
                    dis = abs(bpmRange[idxTmp]-iref)
                    if dis<mnDis:
                        mnDis = dis
                        mnDisIdx = idxTmp
                pk.append([iref, mnDis, bpmRange[mnDisIdx], 20*np.log10(Px[mnDisIdx])])
                #todo handle edge -> assuming argrelmax does not consider edge as pk
            else: # no local pks
                pk.append([iref, 999, 0, 999]) # refHR, dis, HR, Px
            strOut += fmt%(pk[i][0], pk[i][1], pk[i][2], pk[i][3]) + '\n'
        return pk, strOut

    def findMinimumSpectrumPeaks(self):
        '''
        The function uses the spectrum to find the minimum (specMins) and amplitude (specMinVals) of spectrum.
        It finds frequency minimums in the frequency range (match magnitude).
        '''
        maxSpectrum = np.max(self.X)
        self.specMins, self.specMinVals = self.getLocalPeaks(maxSpectrum - self.X, self.dBcutoff,peakLo=self.peakLo_fft,
                                                             peakHi=self.peakHi_fft, searchBpmBins=self.bpmBins_fft,
                                                             numPeaks=self.numPeaks)

        for i in np.arange(np.shape(self.specMins)[0]):
           if self.specMins[i] == 0:
                self.specMinVals[i] = 0
           else:
               self.specMinVals[i] = maxSpectrum - self.specMinVals[i]

        return


if __name__ == '__main__':
    
    
    parser = argparse.ArgumentParser()
    
    #peak filter parameters
    parser.add_argument('--file',help='input file',default=None)
    parser.add_argument('--bpm',help='test signal bpm',type=int,default=60)
    parser.add_argument('--fftLen',help='specify the fft length',type=int,default=8192)
    parser.add_argument('--useWelch',help='use welch''s method to compute power spectrum',
                        default=False, action='store_true')
    args = parser.parse_args()


    ps = PowerSpec(winLen=1000, fftLen=args.fftLen, useWelch=args.useWelch, numHarmonics=3, numPeaks=3)
    f = args.bpm/60.0

    #create a fundmental frequency and the first 2 harmonics
    N = 1000
    x=np.cos(2*np.pi*np.arange(0,N)*f/100) + \
        .5*np.cos(2*np.pi*np.arange(0,N)*2*f/100) + \
        .25*np.cos(2*np.pi*np.arange(0,N)*3*f/100)


    ps.getPowerSpec(x)
    ps.getHarmProdSpec()

    print('Spec Peaks: ', ps.specPeaks)
    print('HPS Peaks : ',ps.hpsFreqPeaks)

    # ps.plotTD()
    # ps.plot()
    pl.show()


