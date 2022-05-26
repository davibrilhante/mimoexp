#! usr/bin/env python3
# -*- coding : utf8 -*-

import numpy as np
from utils import tolin

class Noise:
    def __init__(self, size):
        self.size = size

    def sample(self):
        pass

class Fading:
    def __init__(self, *args, **kwargs):
        pass


    def realization(self):
        pass

class Channel:
    def __init__(self, tx_dim, rx_dim):
        self.tx_dim = tx_dim
        self.rx_dim = rx_dim
        self.noise = None
        self.fading = None

    def setNoise(self, noise):
        pass

    def setFading(self, fading):
        pass

    def response(self, **params):
        pass


class AwgnNoise(Noise):
    def __init__(self, size):
        super().__init__(size)
        self.sigma = 1.0
        self.isComplex = True
    
    def snrdb_to_sigma(self, snr, signal, L=1):
        if not isinstance(signal[0], list):
            sig_power = L*sum(np.absolute(signal)**2)/len(signal)
        else:
            sig_power = L*sum(sum(np.absolute(signal)**2))/len(signal)

        N0 = sig_power/tolin(snr)

        self.sigma = np.sqrt(N0/2)

    def sample(self):
        if self.isComplex:
            nsample = self.sigma*(np.random.standard_normal(self.size) + 
                    1j*np.random.standard_normal(self.size)) #/np.sqrt(2)
        else:
            nsample = self.sigma*np.random.normal(self.size)
            
        return nsample


class RayleighFading(Fading):
    def __init__(self, tx_dim, rx_dim):
        super().__init__(tx_dim, rx_dim)
        self.tx_dim = tx_dim
        self.rx_dim = rx_dim
        self.noise = None

    def realization(self):
        # rayleigh channel sample
        fading_sample = np.random.normal(0,1,(self.rx_dim,self.tx_dim))+1j*(np.random.normal(0,1,(self.rx_dim,self.tx_dim)))

        # Normalization
        for k in range(self.tx_dim):
            fading_sample[:,k] = np.sqrt(1/np.var(fading_sample[:,k]))*fading_sample[:,k]

        return fading_sample

class AwgnRayleighChannel(Channel):
    def __init__(self, tx_dim, rx_dim):
        super().__init__(tx_dim, rx_dim)
        self.noise = AwgnNoise(rx_dim)
        self.fading = RayleighFading(tx_dim, rx_dim)

    def setNoise(self, noise):
        if not issubclass(noise, Noise):
            raise TypeError(': noise is not an object from Noise')

        else:
            self.noise = noise


    def setFading(self, fading):
        if not issubclass(fading, Fading):
            raise TypeError(': fading is not an object from Fading')

        else:
            self.fading = fading

    def response(self, signal, snr, L=1):
        if np.isrealobj(signal):
            self.noise.isComplex = False
        else:
            self.noise.isComplex = True

        channel_sample = []
        noise_var = []


        received = []
        for i in signal:
            channel_sample.append(self.fading.realization())

            signal_at_rx = np.matmul(channel_sample[-1],np.transpose(i))
            self.noise.snrdb_to_sigma(snr, signal_at_rx, L)
            noise_var.append(self.noise.sigma**2)

            received.append(signal_at_rx + self.noise.sample())

        return received, channel_sample, noise_var
