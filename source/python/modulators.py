#! /usr/bin/env python3
# -*- docing : utf8 -*-

import numpy as np
from abc import ABC, abstractmethod

from utils import base_2_to_10_array, binaryToGray
from channel import AwgnRayleighChannel
from demodulators import DemodulatorBPSK
from detectors import ZeroForcing

__all__ = ["Modulator","ModulatorBPSK", "ModulatorQPSK", "ModulatorMPSK"]

class Modulator:
    def __init__(self, order):
        self.mod_order = order
        self.bits_per_symbol = np.log2(order);

    def modulate(self, data_input):
        pass


class ModulatorBPSK(Modulator):
    def __init__(self):
        order = 2
        super().__init__(order)
        #self.mod_order = 2
        #self.bits_per_symbol = 1

    def modulate(self, data_input):
        #Input must always be a tuple of 0 or 1, no matter the type
        modulated = []
        for bit in data_input:
            modulated.append((2*int(bit) - 1) + 0j)

        return modulated

class ModulatorQPSK(Modulator):
    def __init__(self):
        order = 4
        super().__init__(order)
        #self.mod_order = 4
        #self.bits_per_symbol = 2

    def modulate(self, data_input):
        modulated = []
        chunk = []
        for bit in data_input:
            chunk.append(bit)
            if len(chunk) == 2:
                symb = base_2_to_10_array(chunk)
                i_sample = np.cos((np.pi/4)*(2*symb + 1))
                q_sample = np.sin((np.pi/4)*(2*symb + 1))

                modulated.append(i_sample + q_sample*1j)
                chunk = []
        
        return modulated


class ModulatorMPSK(Modulator):
    def __init__(self, order):
        super().__init__(order)
        #self.mod_order = order
        #self.bits_per_symbol = np.log2(order)

    def modulate(self, data_input):
        if self.mod_order == 0:
            raise AttributeError('Modulation order not especified!')

        self.bits_per_symbol = np.log2(self.mod_order)

        modulated = []
        chunk = []
        for bit in data_input:
            chunk += bit
            if len(chunk) == self.bits_per_symbol:
                symb = base_2_to_10_array(chunk)
                i_sample = np.cos((np.pi/self.mod_order)*(symb))
                q_sample = np.sin((np.pi/self.mod_order)*(symb))
                modulated.append(i_sample + q_sample*1j)

                chunk = []
        
        return modulated

class ModulatorMQAM(Modulator):
    '''Based on the python code for QAM modulator found at:
    https://www.gaussianwaves.com/2012/10/qam-modulation-simulation-matlab-python/
    '''
    def __init__(self, order):
        super().__init__(order)
        self.graycode = True

    def modulate(self, data_input):
        if self.mod_order == 0:
            raise AttributeError('Modulation order not especified!')

        self.bits_per_symbol = np.log2(self.mod_order)
        grid_dimension = np.sqrt(self.mod_order).astype(int)

        modulated = []
        chunk = '' 
        for bit in data_input:
            chunk += bit
            if len(chunk) == self.bits_per_symbol:
                if self.graycode:
                    chunk = binaryToGray(chunk)

                symb = base_2_to_10_array(chunk)

                (x, y) = np.divmod(symb, grid_dimension)
                i_sample = 2*x+1-grid_dimension
                q_sample = 2*y+1-grid_dimension

                modulated.append(i_sample + q_sample*1j)

                chunk = ''

        return modulated

    def constellation(self):
        stream = ''

        for num in range(self.mod_order):
            stream += np.binary_repr(num,int(self.bits_per_symbol))

        return self.modulate(stream)

