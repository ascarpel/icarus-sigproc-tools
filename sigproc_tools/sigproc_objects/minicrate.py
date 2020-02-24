import numpy as np
from sigproc_tools.sigproc_functions.noiseAnalysis import getPowerVec
from sigproc_tools.sigproc_functions.noiseProcessing import getPedestalsAndRMS

class MiniCrate:
    """
    Simple minicrate class holding just the portion of the total event.
    Correlations, FFT, etc are first done on the same minicrate
    """

    def __init__ (self, num, waveforms):
        """
        Save the only part of the waveform who belongs to the crate, unsing the
        number of channels per crate and the crate number
        """

        self.nchannels_minicrate = 576
        self.nboard = 64
        self.ncables = 32

        self.ch_per_crate=576
        self._num = num
        self._waveforms = waveforms[num*self.ch_per_crate:(num+1)*self.ch_per_crate]
        self._channels = np.arange(num*self.ch_per_crate, (num+1)*self.ch_per_crate)

        self.plane = self._separatePlanes()

    def _getRMS(self, waveform):
        rms = np.sqrt(np.mean(waveform**2,axis=1))
        return rms

    def _getMedianNoiseCorrection(self, waveforms):
            median = np.median(waveforms,axis=0)
            waveLessCoherent = waveforms - median.transpose()
            return waveLessCoherent

    def _removeCoherentNoise(self, waveforms,grouping):
            # Define placeholders for the output arrays
            waveLessCoherent = np.array([0])

            nChannels = waveforms.shape[0]

            for idx in range(0,nChannels,grouping):
                temp = self._getMedianNoiseCorrection(waveforms[idx:idx+grouping,:])
                if idx == 0:
                    waveLessCoherent = temp
                else:
                    waveLessCoherent = np.concatenate((waveLessCoherent,temp),axis=0)

            return waveLessCoherent


    def _separatePlanes(self):
        """
        Rearrange the waveform to match groups each 32 channels
        """
        groups = np.arange(int(self.nchannels_minicrate/self.ncables))
        plane={
        '0' : np.concatenate([ self._waveforms[self.ncables*n:self.ncables*(n+1)] for n in groups[(groups%2==0)] ]),
        '1' : np.concatenate([ self._waveforms[self.ncables*n:self.ncables*(n+1)] for n in groups[~(groups%2==0)] ])
        }
        return plane

    def remove_coherent_noise(self, grouping):
        self.plane = { key : self._removeCoherentNoise(item, grouping) for key, item in self.plane.items() }

    def get_corr(self):
        corrd = { key : np.corrcoef(item) for key, item in self.plane.items() }
        return corrd

    def get_chanfft(self, maxFrequency):
        fftd = { key : getPowerVec(item, maxFrequency)[1] for key, item in self.plane.items() }
        return fftd

    def get_chanrms(self):
        rmsd = { key : self._getRMS(item) for key, item in self.plane.items() }
        return rmsd
