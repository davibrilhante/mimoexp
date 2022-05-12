#! /usr/bin/env python3
# -*- docing : utf8 -*-

import numpy as np
from abc import ABC, abstractmethod
from channel import AwgnRayleighChannel

from modulators import ModulatorBPSK
from demodulators import DemodulatorBPSK
from detectors import ZeroForcing

if __name__ == '__main__':
    #Byte string to be transmitted
    data = '0001101100011011'

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

    #Signal at the receivers antennas
    received, channel_estimate = channel.response(encoded_data, SNR)

    # Symbol by symbol Zero-Forcing detection
    zf_detector = ZeroForcing()
    constellation = [-1+0j, 1+0j]
    detected = []
    for n, symbol in enumerate(received):
        detected.append(zf_detector.detect(symbol, channel_estimate[n], constellation))

    print('Transmitted data: ',data)
    print('\n\nTransmitted Symbols: ', modulated_data)
    print('\n\nReceived symbols: ', received)
    print('\n\nDetected symbols: ', detected)
    print('\n\nDemodulated data: ', demodulator.demodulate(detected))
