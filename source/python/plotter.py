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
                },
            'zfsic':{
                'file': 'zf-sic',
                'color':'b',
                'label':'ZF-SIC',
                'ls':'-.'
                }
            }


    stream_length = [256, 512, 1024]
    m_order = [2, 4, 16]
    n_antennas = [2,4,8]


    for length in stream_length:
        for n_t in n_antennas:
            for m in m_order:
                fig, ax1 = plt.subplots()
                ax1.set_ylabel('BER', color='r')
                ax1.set_ylim([10e-5,0.75])
                ax1.set_xlabel('SNR [dB]')
                ax1.tick_params(axis='y', colors='r')

                ax2 = ax1.twinx()
                ax2.set_ylim([0.001,0.2])
                ax2.set_ylabel('Time (s)', color='b')
                ax2.tick_params(axis='y', colors='b')

                

                for tec in techniques.keys():
                    filename = 'results/{dir_}/{prefix}-{len_}-{n}-{m}'.format(dir_=tec,prefix=tec,len_=length,n=n_t,m=m) #techniques[tec]['file']
                    print(filename)
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

                plt.title('{m}-QAM {n}x{n} AWGN Rayleigh Channel \n {len_} bit long stream'.format(
                                n=n_t,m=m,len_=length))
                plt.grid()
                ax1.legend()
                fig.tight_layout()
                #plt.show()
                plt.savefig('{m}QAM-{n}x{n}-{len_}'.format(n=n_t,m=m,len_=length))

                plt.figure().clear()
                plt.close()
                plt.cla()
                plt.clf()
