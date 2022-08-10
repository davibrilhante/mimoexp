#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         

import numpy as np                                                              
from abc import ABC, abstractmethod                                             
from utils import grayToBinary

class Demodulator:
    def __init__(self, order):
        self.mod_order = order
        self.bits_per_symbol = int(np.log2(self.mod_order))

    def demodulate(self, data_input):
        pass

class DemodulatorBPSK(Demodulator):
    def __init__(self):
        order = 2
        super().__init__(order)

    def demodulate(self, data_input):
        data_output = ''

        if (isinstance(data_input[0], list) or 
                isinstance(data_input[0], np.ndarray)):
            for stream in data_input:
                for symbol in stream:
                    data_output += str(max(0, int(np.real(symbol)))) 

        else:
            for symbol in data_input:
                data_output += str(max(0, int(np.real(symbol)))) 

        return data_output

class DemodulatorQPSK(Demodulator):
    def __init__(self):
        order = 4
        super().__init__(order)

        self.halfpi = np.pi/2
        self.doublepi = 2*np.pi

    def demodulate(self, data_input):
        data_output = ''

        if (isinstance(data_input[0], list) or 
                isinstance(data_input[0], np.ndarray)):
            for stream in data_input:
                for symbol in stream:
                    angle = np.angle(symbol)
                    if angle < 0:
                        angle += self.doublepi

                    data_output += np.binary_repr(int(np.divmod(angle, self.halfpi)[0]),
                                                    width=self.bits_per_symbol)
        else:
            for symbol in data_input:
                angle = np.angle(symbol)
                if angle < 0:
                    angle += self.doublepi

                data_output += np.binary_repr(int(np.divmod(angle, self.halfpi)[0]),
                                                width=self.bits_per_symbol)
        return data_output


class DemodulatorMPSK(Demodulator):
    def __init__(self, order):
        super().__init__(order)
        self.angle_step = 2*np.pi/self.mod_order

    def demodulate(self, data_input):
        data_output = ''

        if (isinstance(data_input[0], list) or 
                isinstance(data_input[0], np.ndarray)):
            for stream in data_input:
                for symbol in stream:
                    angle = np.angle(symbol)
                    if angle < 0:
                        angle += self.doublepi

                    data_output += np.binary_repr(int(np.divmod(angle, self.angle_step)[0]),
                                                    width=self.bits_per_symbol)
        else:
            for symbol in data_input:
                for symbol in stream:
                    angle = np.angle(symbol)
                    if angle < 0:
                        angle += self.doublepi

                    data_output += np.binary_repr(int(np.divmod(angle, self.angle_step)[0]),
                                                    width=self.bits_per_symbol)

        return data_output


class DemodulatorMQAM(Demodulator):
    def __init__(self, order):
        super().__init__(order)
        self.graycode = True

        self.halfpi = np.pi/2
        self.doublepi = 2*np.pi

    def demodulate(self, data_input):
        data_output = ''
        grid_dimension = np.sqrt(self.mod_order).astype(int)

        if (isinstance(data_input[0], list) or 
                isinstance(data_input[0], np.ndarray)):
            for stream in data_input:
                for symbol in stream:
                    x=(np.real(symbol)+grid_dimension-1)/2
                    y=(np.imag(symbol)+grid_dimension-1)/2

                    data = np.binary_repr(int(x*grid_dimension + y),
                                            width=self.bits_per_symbol)
                    if self.graycode:
                        data = grayToBinary(data)

                    data_output += data
        else:
            for symbol in data_input:
                x=(np.real(symbol)+grid_dimension-1)/2
                y=(np.imag(symbol)+grid_dimension-1)/2

                data = np.binary_repr(int(x*grid_dimension + y),
                                        width=self.bits_per_symbol)

                if self.graycode:
                    data = grayToBinary(data)

                data_output += data

        return data_output
