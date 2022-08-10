#! /usr/bin/env python3
# -*- coding : utf8 -*-

import numpy as np
from matplotlib import pyplot as plt
from ast import literal_eval

if __name__ == '__main__':
    techniques = {
            'mmse':{
                'file': 'linearmmse',
                'color':'r',
                'label':'LMMSE',
                'ls':'--'
                },

            'zf':{
                'file': 'zeroforcing',
                'color':'b',
                'label':'ZF',
                'ls':'-'
                }
            'zf-sic':{
                'file': 'zf-sic',
                'color':'b',
                'label':'ZF-SIC',
                'ls':'-.'
                }
            }

    fig, ax1 = plt.subplots()
    ax1.set_ylabel('BER', color='r')
    ax1.set_ylim([10e-5,0.75])
    ax1.set_xlabel('SNR [dB]')
    ax1.tick_params(axis='y', colors='r')

    ax2 = ax1.twinx()
    ax2.set_ylim([0.01,0.2])
    ax2.set_ylabel('Time (s)', color='b')
    ax2.tick_params(axis='y', colors='b')


    for tec in techniques.keys():
        filename = techniques[tec]['file']
        index = []
        ber_mean = []
        time_mean = []

        with open(filename) as infile:
            for line in infile.readlines():
                stat = literal_eval(line)
                index.append(stat['snr'])
                ber_mean.append(stat['ber_mean'])
                time_mean.append(stat['time_mean'])

        for n, lista in enumerate([ber_mean, time_mean]):
            _index = index.copy()
            zipped = zip(_index, lista)
            sorted_= sorted(zipped)
            tupled = zip(*sorted_)

            _index, lista= [list(Tuple) for Tuple in tupled]


            if n == 0:
                ax1.semilogy(_index, lista, 
                        label=techniques[tec]['label'],
                        color='r',#techniques[tec]['color'],
                        linestyle=techniques[tec]['ls'])
            else:
                ax2.semilogy(_index, lista, 
                        color='b',#techniques[tec]['color'],
                        linestyle=techniques[tec]['ls'])

    plt.title('BPSK 4x4 AWGN Rayleigh Channel') 
    #plt.grid()
    ax1.legend()
    fig.tight_layout()
    plt.show()
