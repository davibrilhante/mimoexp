#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         
import numpy as np                                                              

from abc import ABC, abstractmethod                                             
from channel import AwgnRayleighChannel                                         
from demodulators import DemodulatorBPSK 
from utils import moorePenroseInv
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

    def detect(self, data_input, channel, constellation, noiseVar):

        channel_transf = np.matmul(
                                    np.linalg.pinv(
                                        np.matmul(np.transpose(np.conj(channel)), channel) + 
                                        noiseVar*np.identity(len(data_input))),
                                    np.transpose(np.conj(channel)))

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
        channel_pinv = np.linalg.pinv(channel) #moorePenroseInv(channel) 

        # Line 9c
        channelPower = np.absolute(np.multiply(channel, channel))
        channelPower = np.transpose(channelPower)

        minSnrIndex = min([[n, sum(i)] for n,i in enumerate(channelPower)], key = lambda x : x[1])[0]

        estimate = [0 for i in data_input]
        decoded_symbols = []
        received_symbols = data_input.copy()


        while iteration < len(data_input):
            decoded_symbols.append(minSnrIndex)

            # Line 9d
            minSnrComponent = channel_pinv[:,minSnrIndex]
            
            # Line 9e
            minSnrSymbol = np.matmul(minSnrComponent, received_symbols)
            #print(minSnrSymbol,minSnrComponent,received_symbols)

            # Line 9f
            estimate[minSnrIndex] = min([[i, np.absolute(minSnrSymbol - i)**2] for i in constellation], key = lambda x: x[1])[0]
            #print(minSnrSymbol,estimate[minSnrIndex])

            # Line 9g
            received_symbols = received_symbols - np.multiply(estimate[minSnrIndex],channel[:,minSnrIndex])

            # Line 9h
            channel[:,minSnrIndex] = 0

            channel_pinv = np.linalg.pinv(channel) #moorePenroseInv(channel)

            # Line 9i
            channelPower = np.absolute(np.multiply(channel, channel))
            channelPower = np.transpose(channelPower)
                
            try:
                minSnrIndex = min([[n, sum(i)] for n,i in enumerate(channelPower) if n not in decoded_symbols], key = lambda x : x[1])[0]
            except ValueError:
                break

            # Line 9j
            iteration += 1

        return estimate
