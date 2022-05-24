#! /usr/bin/env python3
# -*- docing : utf8 -*-

import numpy as np
from abc import ABC, abstractmethod
from channel import AwgnRayleighChannel

from utils import gen_binary_str, bit_error
from modulators import ModulatorBPSK
from demodulators import DemodulatorBPSK
from detectors import ZeroForcing, LinearMse, ZeroForcingSic

if __name__ == '__main__':
    #Byte string to be transmitted
    data_len = 1024
    data = gen_binary_str(data_len)

    #MIMO device specificantions
    tx_dimension = 2
    rx_dimension = 2

    #desired SNR, must vary on a loop
    SNR = 10

    #BPSK modulation example
    modulator = ModulatorBPSK()
    demodulator = DemodulatorBPSK()
    modulated_data = modulator.modulate(data)


    #Reshape the signal to the number of transmitters antennas
    encoded_data = np.reshape(modulated_data,(-1,rx_dimension))


    #Creating the channel to add interference and noise
    channel = AwgnRayleighChannel(tx_dimension,rx_dimension)
    
    #Set the noise configuration
    channel.noise.snrdb_to_sigma(SNR,modulated_data)
    noisevar = channel.noise.sigma**2

    #Signal at the receivers antennas
    received, channel_estimate = channel.response(encoded_data, SNR)

    # Symbol by symbol Zero-Forcing detection
    zf_detector = ZeroForcing()
    constellation = [-1+0j, 1+0j]
    zf_detected = []
    for n, symbol in enumerate(received):
        zf_detected.append(zf_detector.detect(symbol, channel_estimate[n], constellation))

    zf_demodulated =  demodulator.demodulate(zf_detected)

    lmse_detector = LinearMse()
    lmse_detected = []
    for n, symbol in enumerate(received):
        lmse_detected.append(lmse_detector.detect(
                        symbol, channel_estimate[n], constellation, noisevar))

    lmse_demodulated = demodulator.demodulate(lmse_detected)
    
    # Symbol by symbol Zero-Forcing detection
    sic_detector = ZeroForcingSic()
    constellation = [-1+0j, 1+0j]
    sic_detected = []
    for n, symbol in enumerate(received):
        sic_detected.append(sic_detector.detect(symbol, channel_estimate[n], constellation))

    sic_demodulated =  demodulator.demodulate(sic_detected)

    print('Transmitted data: ',data)
    #print('\n\nTransmitted Symbols: ', modulated_data)
    #print('\n\nReceived symbols: ', received)
    #print('\n\nDetected symbols: ', detected)
    print('\n\nDemodulated data: ', zf_demodulated)
    print('\n\nZF Bit Error: ',bit_error(data, zf_demodulated))
    print('\nLMSE Bit Error: ',bit_error(data, lmse_demodulated))
    print('\nZF-SIC Bit Error: ',bit_error(data, sic_demodulated))
