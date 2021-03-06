#! /usr/bin/env python3
# -*- docing : utf8 -*-

import numpy as np
from abc import ABC, abstractmethod

__all__ = ['base_2_to_10_array', 'gen_binary_str', 'gen_binary_array', 'bit_error',
        'bit_error_rate']

def base_2_to_10_array(arr):#, i1, i2):
    res = 0
    for bit in arr:#[i1:i2][::-1]:
        bit = int(bit)
        res = (res << 1) ^ bit

    return res

def gen_binary_str(length : int)->str:
    return ''.join([np.random.choice(['0','1']) for i in range(length)])

def gen_binary_array(length : int)->list:
    return [np.random.choice(['0','1']) for i in range(length)]

def bit_error(orig, recv)->int:
    error = 0
    for i, j in zip(orig, recv):
        if i != j:
            error += 1

    return error

def bit_error_rate(orig, recv)->float:
    error = bit_error(orig, recv)
    return error/len(orig)


def tolin(s):
    return 10**(0.1*s)

def binaryToGray(b):
    b = int(b,2)
    b ^= (b>>1)

    return bin(b)[2:]

def grayToBinary(b):
    """Convert Gray codeword to binary and return it."""
    b = int(b, 2) # convert to int
 
    mask = b
    while mask != 0:
        mask >>= 1
        b ^= mask
 
    return bin(b)[2:]
