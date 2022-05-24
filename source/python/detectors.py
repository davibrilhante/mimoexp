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
        channel_pinv = moorePenroseInv(channel)

        tildey = np.dot(channel_pinv, np.transpose(data_input))

        estimate = []
        for symbol in tildey:
            approx = min([[i, np.absolute(symbol - i)**2] for i in constellation], key= lambda x: x[1])
            estimate.append(approx[0])

        return estimate

class LinearMse(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation, noiseVar):
        channel_transf = np.dot(
                np.linalg.inv(
                        np.dot(np.conj(channel), channel)+noiseVar*np.identity(len(data_input))
                    ),np.conj(channel))

        tildey = np.dot(channel_transf, np.transpose(data_input))

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

        print(data_input)

        # Line 9a, but starting from 0
        iteration = 0

        # Line 9b
        channel_pinv = moorePenroseInv(channel) 

        # Line 9c
        channelPower = []
        for i in channel_pinv:
            channelPower.append(np.absolute(i)**2)
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
            minSnrSymbol = np.dot(minSnrComponent, received_symbols)

            # Line 9f
            estimate[minSnrIndex] = min([[i, np.absolute(minSnrSymbol - i)**2] for i in constellation], key = lambda x: x[1])[0]
            #print(minSnrSymbol,estimate[minSnrIndex])

            # Line 9g, operated over every column except the already decoded
            for n, symbol in enumerate(data_input):
                if n not in decoded_symbols: 
                    received_symbols[n] = received_symbols[n] - estimate[minSnrIndex]*channel[n,minSnrIndex]

            # Line 9h
            #channel_pinv[:,minSnrIndex] = 0
            channel[:,minSnrIndex] = 0

            channel_pinv = moorePenroseInv(channel)

            # Line 9i
            channelPower = np.absolute(channel_pinv)**2
            for n in decoded_symbols:
                # Will eliminate the already decoded symbol columns for the argmin
                channelPower[:,n] = float('inf')
            channelPower = np.transpose(channelPower)
                
            minSnrIndex = min([[n, sum(i)] for n,i in enumerate(channelPower)], key = lambda x : x[1])[0]

            # Line 9j
            iteration += 1

        return estimate
