#! /usr/bin/env python3
# -*- docing : utf8 -*-

import numpy as np

from channel import AwgnRayleighChannel

from utils import gen_binary_str, bit_error, bit_error_rate
from modulators import ModulatorMQAM
from demodulators import DemodulatorMQAM
from detectors import ZeroForcing
from detectors import LinearMmse
from detectors import ZeroForcingSic
from detectors import MeanSquaredErrorSic
from detectors import SphereDetector 
from precoders import LatticeReductionAided

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
parser.add_argument('-m','--morder',type=int,default=4,help='QAM Modulation Order')


group = parser.add_mutually_exclusive_group()
group.add_argument('--zf',action='store_true')
group.add_argument('--mmse',action='store_true')
group.add_argument('--zfsic',action='store_true')
group.add_argument('--lra',action='store_true')
group.add_argument('--sd',action='store_true')
group.add_argument('--msesic',action='store_true')

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

    #high.start_counters([events.PAPI_FP_OPS,])
    #high.flops()

    for i in range(args.executions):
        start = default_timer()

        np.random.seed(i)

        data = gen_binary_str(data_len)

        #MQAM modulation example
        modulator = ModulatorMQAM(args.morder)
        modulator.graycode = False
        demodulator = DemodulatorMQAM(args.morder)
        demodulator.graycode = False
        modulated_data = modulator.modulate(data)

        #Reshape the signal to the number of transmitters antennas
        encoded_data = np.reshape(modulated_data,(-1,tx_dimension))


        #Creating the channel to add interference and noise
        channel = AwgnRayleighChannel(tx_dimension,rx_dimension)
        
        #Signal at the receiver antennas
        received, channel_estimate, noisevar = channel.response(encoded_data, SNR)
        sqrth = np.sqrt(2)/2
        #constellation = [-1-1j, -1+1j, 1-1j,1+1j] #np.multiply(sqrth,[1+1j,-1+1j,-1-1j,1-1j]) #[-1+0j, 1+0j]
        constellation = modulator.constellation()

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
                #print(noisevar[n], lmse_detector.snrdb_to_sigma(SNR,symbol))
                lmse_detected.append(lmse_detector.detect(
                                symbol, channel_estimate[n], constellation, SNR))
                                #symbol, channel_estimate[n], constellation, noisevar[n]))

            demodulated = demodulator.demodulate(lmse_detected)

        elif args.zfsic:    
            # Symbol by symbol Zero-Forcing detection
            sic_detector = ZeroForcingSic()
            sic_detected = []
            for n, symbol in enumerate(received):
                sic_detected.append(sic_detector.detect(symbol, channel_estimate[n], constellation))

            demodulated =  demodulator.demodulate(sic_detected)

        elif args.lra:    
            precoder = LatticeReductionAided()
            precoder.type = 'complex'
            received, channel_estimate, noisevar = channel.response(encoded_data, SNR)
            U, h_red = precoder.precode(channel_estimate)

            for i in range(len(encoded_data)):
                signal = np.matmul(np.matmul(channel_estimate[i],U[i]),encoded_data[i])
                channel.noise.snrdb_to_sigma(SNR, [signal])
                received[i] =  signal + channel.noise.sample()

            # Symbol by symbol Zero-Forcing detection
            detector = ZeroForcing()
            detected = []
            for n, symbol in enumerate(received):
                detected.append(detector.detect(symbol, h_red[n], constellation))

            demodulated =  demodulator.demodulate(detected)

        elif args.msesic:
            lmse_detector = MeanSquaredErrorSic()
            lmse_detected = []
            for n, symbol in enumerate(received):
                lmse_detected.append(lmse_detector.detect(
                                symbol, channel_estimate[n], constellation, SNR))
                                #symbol, channel_estimate[n], constellation, noisevar[n]))

            demodulated = demodulator.demodulate(lmse_detected)

        elif args.sd:    
            # Symbol by symbol Sphere Detector detection
            sd_detector = SphereDetector()
            sd_detected = []
            for n, symbol in enumerate(received):
                sd_detected.append(sd_detector.detect(symbol, channel_estimate[n], constellation, SNR))

            demodulated =  demodulator.demodulate(sd_detected)

        stop = default_timer()

        output['error'].append(bit_error(data, demodulated))
        output['ber'].append(bit_error_rate(data, demodulated))
        output['time'].append(stop - start)

    #nflops=high.flops()
    #high.stop_counters()

    final = {}
    final['snr'] = args.snr
    final['ber_mean'] = np.mean(output['ber'])
    final['ber_std'] = np.std(output['ber'])
    final['time_mean'] = np.mean(output['time'])
    final['time_std'] = np.std(output['time'])
    #final['flops'] = nflops

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
