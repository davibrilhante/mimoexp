#! /usr/bin/env python3
# -*- docing : utf8 -*-

import numpy as np
from channel import AwgnRayleighChannel

from utils import gen_binary_str, bit_error, bit_error_rate
from modulators import ModulatorQPSK
from demodulators import DemodulatorQPSK
from detectors import ZeroForcing, LinearMmse, ZeroForcingSic

from argparse import ArgumentParser
from timeit import default_timer
from json import dumps

parser = ArgumentParser()
parser.add_argument('-l','--length',type=int,default=1024,help='Random binary str input data length')
parser.add_argument('-t','--tx',type=int,default=2,help='Tx channel length')
parser.add_argument('-r','--rx',type=int,default=2,help='Rx channel length')
parser.add_argument('-S','--snr',type=int,default=10,help='Signal SNR')
parser.add_argument('-s','--seed',type=int,default=1,help='random seed')
parser.add_argument('-x','--executions',type=int,default=1000,help='how many realizations')


group = parser.add_mutually_exclusive_group()
group.add_argument('--zf',action='store_true')
group.add_argument('--mmse',action='store_true')
group.add_argument('--zfsic',action='store_true')

args = parser.parse_args()



if __name__ == '__main__':
    #np.random.seed(args.seed)


    #Byte string to be transmitted
    data_len = args.length

    #MIMO device specificantions
    tx_dimension = args.tx
    rx_dimension = args.rx


    #desired SNR, must vary on a loop
    SNR = args.snr

    output = {}
    output['error'] = []
    output['ber'] = []
    output['time'] = []

    for i in range(args.executions):
        start = default_timer()

        np.random.seed(i)

        data = gen_binary_str(data_len)

        #BPSK modulation example
        modulator = ModulatorQPSK()
        demodulator = DemodulatorQPSK()
        modulated_data = modulator.modulate(data)


        #Reshape the signal to the number of transmitters antennas
        encoded_data = np.reshape(modulated_data,(-1,tx_dimension))


        #Creating the channel to add interference and noise
        channel = AwgnRayleighChannel(tx_dimension,rx_dimension)
        
        #Signal at the receivers antennas
        received, channel_estimate, noisevar = channel.response(encoded_data, SNR)
        sqrth = np.sqrt(2)/2
        constellation = np.multiply(sqrth,[1+1j,-1+1j,-1-1j,1-1j]) #[-1+0j, 1+0j]

        # Symbol by symbol Zero-Forcing detection
        if args.zf:
            zf_detector = ZeroForcing()
            zf_detected = []
            for n, symbol in enumerate(received):
                zf_detected.append(zf_detector.detect(symbol, channel_estimate[n], constellation))

            demodulated =  demodulator.demodulate(zf_detected)

        elif args.mmse:
            lmse_detector = LinearMmse()
            lmse_detected = []
            for n, symbol in enumerate(received):
                lmse_detected.append(lmse_detector.detect(
                                symbol, channel_estimate[n], constellation, noisevar[n]))

            demodulated = demodulator.demodulate(lmse_detected)

        elif args.zfsic:    
            # Symbol by symbol Zero-Forcing detection
            sic_detector = ZeroForcingSic()
            sic_detected = []
            for n, symbol in enumerate(received):
                sic_detected.append(sic_detector.detect(symbol, channel_estimate[n], constellation))

            demodulated =  demodulator.demodulate(sic_detected)

        stop = default_timer()

        output['error'].append(bit_error(data, demodulated))
        output['ber'].append(bit_error_rate(data, demodulated))
        output['time'].append(stop - start)

    final = {}
    final['snr'] = args.snr
    final['ber_mean'] = np.mean(output['ber'])
    final['ber_std'] = np.std(output['ber'])
    final['time_mean'] = np.mean(output['time'])
    final['time_std'] = np.std(output['time'])

    y = dumps(final)
    print(y)



    #print(bit_error(data, demodulated))
    #print(bit_error_rate(data, demodulated))
    #print(stop - start)

    #print('Transmitted data: ',data)
    #print('\n\nTransmitted Symbols: ', modulated_data)
    #print('\n\nReceived symbols: ', received)
    #print('\n\nDetected symbols: ', detected)
    #print('\n\nDemodulated data: ', zf_demodulated)
    #print('\n\nZF Bit Error: ',bit_error(data, zf_demodulated))
    #print('\nLMSE Bit Error: ',bit_error(data, lmse_demodulated))
    #print('\nZF-SIC Bit Error: ',bit_error(data, sic_demodulated))
