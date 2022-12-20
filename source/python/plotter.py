#! /usr/bin/env python3
# -*- coding : utf8 -*-

import numpy as np
from matplotlib import pyplot as plt
from ast import literal_eval

plots = {
        'ber-snr-256':{
            'data':'ber',
            'xaxis':'snr',
            'nant':2,
            'mord':[2,4,16],
            'streamlen':256
            },
        'ber-snr-512':{
            'data':'ber',
            'xaxis':'snr',
            'nant':2,
            'mord':[2,4,16],
            'streamlen':512
            },
        'ber-snr-1024':{
            'data':'ber',
            'xaxis':'snr',
            'nant':2,
            'mord':[2,4,16],
            'streamlen':1024
            },
        'time-nant':{
            'xlabel':'Number of Antennas',
            'data':'time',
            'xaxis':'nant',
            'nant':[2,4,8],
            'mord':2,
            'streamlen':1024
            },
        'time-streamlen':{
            'xlabel':'Stream Length [Bits]',
            'data':'time',
            'xaxis':'streamlen',
            'nant':2,
            'mord':2,
            'streamlen':[256, 512,1024]
            },
        'time-mord':{
            'xlabel':'Modulation Order',
            'data':'time',
            'xaxis':'mord',
            'nant':2,
            'mord':[2,4,16],
            'streamlen':1024
            }
        }

techniques = {
        'zf':{
            'file': 'zeroforcing',
            'color':'b',
            'label':'ZF',
            'ls':'-'
            },
        'mmse':{
            'file': 'linearmmse',
            'color':'r',
            'label':'LMMSE',
            'ls':'--'
            },
        'zfsic':{
            'file': 'zf-sic',
            'color':'b',
            'label':'ZF-SIC',
            'ls':'-.'
            },
        'lra':{
            'file': 'lra',
            'color':'g',
            'label':'LRA-SIC',
            'ls':':'
            }
        }

modulation = {
        2:{
            'color':'b',
            'label':'BPSK'
            },
        4:{
            'color':'g',
            'label':'4QAM'
            },
        16:{
            'color':'r',
            'label':'16QAM'
            }
        }


if __name__ == '__main__':
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=13)


    stream_length = [256, 512, 1024]
    m_order = [2, 4, 16]
    n_antennas = [2,4,8]


    data = {}
    for tec in techniques.keys():
        data[tec] = {'ber':{},'time':{}}
        for length in stream_length:
            for n_t in n_antennas:
                for m in m_order:
                    filename = 'results/{dir_}/{prefix}-{len_}-{n}-{m}'.format(dir_=tec,prefix=tec,len_=length,n=n_t,m=m) #techniques[tec]['file']

                    data[tec]['ber'][(length,n_t,m)] = {}
                    data[tec]['time'][(length,n_t,m)] = {}

                    with open(filename) as infile:
                        for line in infile.readlines():
                            stat = literal_eval(line)
                            for type_ in ['ber','time']:
                                data[tec][type_][(length,n_t,m)][stat['snr']] = stat[type_+'_mean']



    for plot in plots.values():
        fig, ax = plt.subplots()
        if plot['xaxis'] == 'snr':
            ax.set_ylabel('BER')
            ax.set_ylim([10e-5,0.75])
            ax.set_xlabel('SNR [dB]')

            for m in plot['mord']:
                for tec in techniques.keys():
                    _index= list(data[tec]['ber'][(plot['streamlen'],plot['nant'],m)].keys())
                    lista = list(data[tec]['ber'][(plot['streamlen'],plot['nant'],m)].values())

                    listline = ax.semilogy(_index, lista, 
                            label=techniques[tec]['label']+', '+modulation[m]['label'],
                            color=modulation[m]['color'],
                            linewidth=0.8,
                            linestyle=techniques[tec]['ls'])

                    if tec == 'zfsic':
                        listline[0].set_dashes((15,10,15,10))
                    if tec == 'mmse':
                        listline[0].set_dashes((15,3,3,3))

            plt.title('{n}x{n} AWGN Rayleigh Channel \n {len_} bit long stream'.format(
                            n=plot['nant'],len_=plot['streamlen']))
            plt.grid()
            ax.legend()
            fig.tight_layout()
            plt.savefig('ber-snr-{n}x{n}-{len_}'.format(n=plot['nant'],m=m,len_=plot['streamlen']))
            
        else:

            for tec in techniques.keys():
                _index = []
                lista = []
                for x in plot[plot['xaxis']]:
                    if plot['xaxis'] == 'nant':
                        length = plot['streamlen']
                        n_t = x
                        m = plot['mord']
                    elif plot['xaxis'] == 'mord':
                        length = plot['streamlen']
                        n_t = plot['nant']
                        m = x
                    elif plot['xaxis'] == 'streamlen':
                        length = x
                        n_t = plot['nant']
                        m = plot['mord']

                    _index.append(x)
                    print(tec, length, n_t, m)
                    try:
                        lista.append(data[tec]['time'][(length,n_t,m)][10])
                    except KeyError:
                        continue
                    #lista = list(data[tec]['time'][(length,n_t,m)].values())

                    if plot['xaxis'] == 'mord':
                        listline = ax.semilogx(_index, lista, 
                                label=techniques[tec]['label'] if x==2 else '',#+', '+modulation[m]['label'],
                                #color=modulation[m]['color'],
                                color='b',
                                linewidth=0.8,
                                linestyle=techniques[tec]['ls'], base=2)
                    else:
                        listline = ax.plot(_index, lista, 
                                label=techniques[tec]['label'] if (x==2 or x==256) else '',#+', '+modulation[m]['label'],
                                #color=modulation[m]['color'],
                                color='b',
                                linewidth=0.8,
                                linestyle=techniques[tec]['ls'])

                    if tec == 'zfsic':
                        listline[0].set_dashes((15,10,15,10))
                    if tec == 'mmse':
                        listline[0].set_dashes((15,3,3,3))

            ax.set_ylabel('Execution time')
            ax.set_xlabel(plot['xlabel'])
            plt.title('AWGN Rayleigh Channel')
            plt.grid()
            ax.legend()
            fig.tight_layout()
            plt.savefig('time-{metric}'.format(metric=plot['xaxis']))




        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()
