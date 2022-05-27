#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         
import numpy as np                                                              

from abc import ABC, abstractmethod                                             
from channel import AwgnRayleighChannel                                         
from demodulators import DemodulatorBPSK 
#from utils import moorePenroseInv
from utils import tolin
#from modulators import *

class Detector:
    def __init__(self):
        pass

    def detect(self, data_input, *args, **kwargs):
        pass

class ZeroForcing(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation):
        channel_pinv = np.linalg.pinv(channel) #moorePenroseInv(channel)

        tildey = np.matmul(channel_pinv, np.transpose(data_input))

        estimate = []
        for symbol in tildey:
            approx = min([[i, np.absolute(symbol - i)**2] for i in constellation], key= lambda x: x[1])
            estimate.append(approx[0])

        return estimate

class LinearMmse(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation, snr): #noiseVar):
        f1 = np.matmul(np.transpose(np.conj(channel)), channel)
        #f2 = 2*noiseVar*np.identity(channel.shape[1])
        f2 = (1/tolin(snr))*np.identity(channel.shape[1])

        channel_transf = np.matmul(np.linalg.pinv(f1+f2),np.transpose(np.conj(channel)))

        tildey = np.matmul(channel_transf, np.transpose(data_input))

        estimate = []
        for symbol in tildey:
            approx = min([[i, np.absolute(symbol - i)**2] for i in constellation], key= lambda x: x[1])
            estimate.append(approx[0])

        return estimate



class ZeroForcingSic(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation):
        # Based on:  
        # V-BLAST: An Architecture for Realizing Very High Data Rates
        #       Over the Rich-Scattering Wireless Channel
        #
        # P. W. Wolniansky, G. J. Foschini, G. D. Golden, R. A. Valenzuela

        # Line 9a, but starting from 0
        iteration = 0


        # Line 9b
        channel_pinv = np.linalg.pinv(channel).T #moorePenroseInv(channel) 
        todecode = [n for n in range(channel.shape[1])]

        # Line 9c
        minSnrIndex = min([[n, np.linalg.norm(channel_pinv[:,n])] for n in todecode], key = lambda x : x[1])[0]


        estimate = [0 for i in data_input]
        decoded_symbols = []
        received_symbols = data_input.copy()


        while iteration < channel.shape[1]:
            # Line 9d
            minSnrComponent = channel_pinv[:,minSnrIndex]
            
            # Line 9e
            minSnrSymbol = np.matmul(minSnrComponent.T, received_symbols)

            # Line 9f
            estimate[minSnrIndex] = min([[i, np.absolute(minSnrSymbol - i)**2] for i in constellation], key = lambda x: x[1])[0]
            decoded_symbols.append(minSnrIndex)

            todecode.remove(minSnrIndex)

            # Line 9g
            received_symbols = received_symbols - np.multiply(channel[:,minSnrIndex],estimate[minSnrIndex])

            # Line 9h
            channel[:,minSnrIndex] = 0

            channel_pinv = np.linalg.pinv(channel).T #moorePenroseInv(channel)

            # Line 9i
            try:
                minSnrIndex = min([[n, np.linalg.norm(channel_pinv[:,n])] for n in todecode], key = lambda x : x[1])[0]
            except ValueError:
                break
            # Line 9j
            iteration += 1

        return estimate
