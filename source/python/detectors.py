#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         
import numpy as np                                                              

from abc import ABC, abstractmethod                                             
from channel import AwgnRayleighChannel                                         
from demodulators import DemodulatorBPSK 
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
        channel_pinv = np.dot(np.linalg.inv(np.dot(np.conj(channel), channel)), np.conj(channel))

        tildey = np.dot(channel_pinv, np.transpose(data_input))

        estimate = []
        for symbol in tildey:
            approx = min([[i, np.absolute(symbol - i)**2] for i in constellation], key= lambda x: x[1])
            estimate.append(approx[0])

        return estimate
