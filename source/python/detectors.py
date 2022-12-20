#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         
import numpy as np                                                              

from abc import ABC, abstractmethod                                             
from channel import AwgnRayleighChannel                                         
from demodulators import DemodulatorBPSK 
#from utils import moorePenroseInv
from utils import tolin
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
        channel_pinv = np.linalg.pinv(channel) #moorePenroseInv(channel)

        tildey = np.matmul(channel_pinv, np.transpose(data_input))

        estimate = []
        for symbol in tildey:
            approx = min([[i, np.absolute(symbol - i)**2] for i in constellation], key= lambda x: x[1])
            estimate.append(approx[0])

        return estimate

class LinearMmse(Detector):
    def __init__(self):
        super().__init__()

    def snrdb_to_sigma(self, snr, data_input, L=1):
        if not isinstance(data_input[0], list):
            sig_power = L*sum(np.absolute(data_input)**2)/len(data_input)
        else:
            sig_power = L*sum(sum(np.absolute(data_input)**2))/len(data_input)

        #print(sig_power, tolin(snr))

        N0 = sig_power/(tolin(snr+1))

        return np.sqrt(N0/2)

    def detect(self, data_input, channel, constellation, snr): #noiseVar):
        f1 = np.matmul(np.transpose(np.conj(channel)), channel)

        #f2 = 2*noiseVar*np.identity(channel.shape[1])
        f2 = (channel.shape[1]/tolin(snr))*np.identity(channel.shape[1])
        #f2 = (2*channel.shape[0]/channel.shape[1]*tolin(-1*snr))*np.identity(channel.shape[1])
        #mod_order = np.log2(len(constellation))
        #f2 = self.snrdb_to_sigma(snr, data_input, mod_order)*np.identity(channel.shape[1])

        channel_transf = np.matmul(np.linalg.pinv(f1+f2),np.transpose(np.conj(channel)))

        tildey = np.matmul(channel_transf, np.transpose(data_input))

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

        # Line 9a, but starting from 0
        iteration = 0


        # Line 9b
        channel_pinv = np.linalg.pinv(channel).T #moorePenroseInv(channel) 
        todecode = [n for n in range(channel.shape[1])]

        # Line 9c
        minSnrIndex = min([[n, np.linalg.norm(channel_pinv[:,n])] for n in todecode], key = lambda x : x[1])[0]


        estimate = [0 for i in data_input]
        decoded_symbols = []
        received_symbols = data_input.copy()


        while iteration < channel.shape[1]:
            # Line 9d
            minSnrComponent = channel_pinv[:,minSnrIndex]
            
            # Line 9e
            minSnrSymbol = np.matmul(minSnrComponent.T, received_symbols)

            # Line 9f
            estimate[minSnrIndex] = min([[i, np.absolute(minSnrSymbol - i)**2] for i in constellation], key = lambda x: x[1])[0]
            decoded_symbols.append(minSnrIndex)

            todecode.remove(minSnrIndex)

            # Line 9g
            received_symbols = received_symbols - np.multiply(channel[:,minSnrIndex],estimate[minSnrIndex])

            # Line 9h
            channel[:,minSnrIndex] = 0

            channel_pinv = np.linalg.pinv(channel).T #moorePenroseInv(channel)

            # Line 9i
            try:
                minSnrIndex = min([[n, np.linalg.norm(channel_pinv[:,n])] for n in todecode], key = lambda x : x[1])[0]
            except ValueError:
                break
            # Line 9j
            iteration += 1

        return estimate


class MeanSquaredErrorSic(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation, snr):
        # Based on:  
        # MMSE Extension of V-BLAST based on Sorted QR Decomposition 
        #
        # Dirk Wubben, Ronald Bohnke, Volker Kuhn, and Karl-Dirk Kammeyer 

        estimate = [0 for i in data_input]

        # Table 1 - line 1 
        CHANNEL = channel
        
        f2 = (channel.shape[1]/tolin(snr))*np.identity(channel.shape[1])
        CHANNEL = np.concatenate((CHANNEL,f2))
        #for i in f2:
        #    CHANNEL.append(i)

        R = np.zeros(channel.shape) 
        #[[0 for j in range(channel.shape[1])] for i in range(channel.shape[0]+channel.shape[1])]


        todecode = [n for n in range(channel.shape[1])]
        norm = []


        # Table 1 - lines 2-4
        for i in range(channel.shape[1]):
            norm.append(np.linalg.norm(CHANNEL[:,i])**2)

        
        # Table 1 - lines 5, 6
        for i in range(channel.shape[1]):
            k_i = min(zip(todecode, norm), key = lambda x : x[1])[0]

            # Table 1 - line 7
            tmp = todecode[i]
            todecode[i] = todecode[k_i]
            todecode[k_i] = tmp

            tmp = norm[i]
            norm[i] = norm[k_i]
            norm[k_i] = tmp

            tmp = R[:,i]
            R[:,i] = R[:,k_i]
            R[:,k_i] = tmp

            limit = channel.shape[0] + i - 1

            for j in range(limit):
                tmp = CHANNEL[j][i]
                CHANNEL[j][i] = CHANNEL[j][k_i]
                CHANNEL[j][k_i] = tmp

            # Table 1 - lines 8,9
            r_i_i = np.sqrt(norm[i])
            CHANNEL[:,i] /= r_i_i

            # Table 1 - lines 10 - 15
            for k in range(i+1,channel.shape[1]):
                r_i_k = np.matmul(np.transpose(np.conj(CHANNEL[:,i])), CHANNEL[:,k])
                CHANNEL[:,k]=CHANNEL[:,k] - r_i_k*CHANNEL[:,i]
                norm[k]=norm[k] - (r_i_k**2)


        k_min = channel.shape[1]

        Q1 = CHANNEL[:channel.shape[0],:]
        Q2 = CHANNEL[channel.shape[0]:,:]
        error = []

        for i in reversed(range(1,channel.shape[1])):
            for l in range(i):
                error.append(np.linalg.norm(Q2[l,:i])**2)

            k_ii = min(zip(range(i),error), key= lambda x : x[1] )[0]

            k_min = min(k_min, k_ii)

            if k_ii < i:
                tmp = Q2[:,i]
                Q2[:,i] = Q2[:,k_ii]
                Q2[:,k_ii] = tmp

                tmp = todecode[i]
                todecode[i] = todecode[k_ii]
                todecode[k_ii] = todecode[i]

            if k_min < i:
                a = Q2[i,k_min:i]
                e = np.concatenate((np.zeros(i-k_min),[1]))

                u = (a - np.linalg.norm(a)*e)/np.linalg.norm(a - np.linalg.norm(a)*e)

                w1 = np.matmul(u,np.transpose(np.conj(a)))
                w2 = np.matmul(a,np.transpose(np.conj(u)))
                w = w1/w2

                householder=np.eye(i-k_min)-(1+w)*np.matmul(np.transpose(np.conj(u)),u)

                Q2[:i,kmin:i+1] = np.matmul(Q2[:i,kmin:i+1], householder)
                Q1[:,kmin:i+1] = np.matmul(Q1[:,kmin:i+1],householder)

        #R = 1/(channel.shape[1]/tolin(snr))*np.linalg.pinv(Q2)

        estimate = np.matmul(np.transpose(np.conj(Q1)),data_input)




        return estimate


class SphereDetector(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation, snr):
        rx_dimension, tx_dimension = np.shape(channel)

        #eigen_values, eigen_vectors = np.linalg.eig(channel)
        eigen_values = np.linalg.eigvals(channel)

        #Initial radius
        #radius = min(eigen_values)
        radius = min([np.linalg.norm(channel[:,i]) for i in range(tx_dimension)])
 
        #QR matrix decomposition
        matrixQ, matrixR = np.linalg.qr(channel)

        #Q matrix split
        matrixQ1 = matrixQ[:,0:tx_dimension]
        matrixQ2 = matrixQ[:,tx_dimension:]

        #calculate s_hat, the estimate
        #estimate = []
        #rough = np.matmul(np.linalg.pinv(channel),data_input)
        #for symbol in rough:
        #    estimate.append(min([[i, np.absolute(symbol - i)**2] for i in constellation], key= lambda x: x[1])[0])
        estimate = np.matmul(np.linalg.pinv(channel),data_input)

        new_radius = [0 for i in range(tx_dimension)]
        new_radius[tx_dimension-1] = np.sqrt(
                radius**2-np.linalg.norm(data_input)**2
                +np.linalg.norm(np.matmul(channel,estimate))**2)

        rel_estimate = [0 for i in range(tx_dimension+1)]
        rel_estimate[tx_dimension] = estimate[tx_dimension-1]

        k = tx_dimension-1
        solution = [0 for i in range(tx_dimension)]

        #step = np.exp(1j*2*np.pi*theta/len(constellation))
        step = 0

        boundflag = True
        while 1:
            if boundflag:
                #step 2
                z = new_radius[k]/matrixR[k][k]
                l_bound = np.ceil(np.real(-z+rel_estimate[k+1])) + 1j*np.ceil(np.imag(-z+rel_estimate[k+1]))
                u_bound = np.floor(np.real(z+rel_estimate[k+1])) + 1j*np.floor(np.imag(z+rel_estimate[k+1]))
                print(new_radius[k], l_bound, u_bound)

                solution[k] = l_bound - step
                

            #step 3
            solution[k] = solution[k]+step

            if solution[k] <= u_bound:
                #step 5
                if k == 0:
                    #step 6
                    boundflag = False

                else:
                    k = k - 1
                    rel_estimate[k+1] = estimate[k] - sum( 
                            (matrixR[k][j]/matrixR[k][k])*(solution[j] - estimate[j]) 
                            for j in range(k+1,tx_dimension-1))

                    new_radius[k] = np.sqrt(
                            new_radius[k+1]**2 -
                            matrixR[k+1][k+1]**2*(solution[k+1] - rel_estimate[k+2])
                            )
                    boundflag = True
                    

            else:
                #step 4
                k = k+1
                if k==tx_dimension:
                    return solution
                else:
                    boundflag = False


class SphereDetector2(Detector):
    def __init__(self):
        super().__init__()

    def detect(self, data_input, channel, constellation, snr):
        rx_dimension, tx_dimension = np.shape(channel)

        #eigen_values, eigen_vectors = np.linalg.eig(channel)
        eigen_values = np.linalg.eigvals(channel)

        #Initial radius
        #radius = tx_dimension*(1/snr)
        radius = min(eigen_values)

        #QR matrix decomposition
        matrixQ, matrixR = np.linalg.qr(channel)

        #Q matrix split
        matrixQ1 = matrixQ[:,0:tx_dimension]
        matrixQ2 = matrixQ[:,tx_dimension:]

        #Least-squares transmitted symbol estimative
        tx_symbol_estimate = np.matmul(np.linalg.pinv(channel), data_input)
        

        iteration_estimate = tx_symbol_estimate[tx_dimension - 1]
        i = tx_dimension-1

        #
        y = np.matmul(np.transpose(np.conj(matrixQ1)), data_input)
        iter_y = y[i]

        #iter_radius = np.sqrt(radius**2 - 
        #        np.linalg.norm(np.matmul(np.transpose(np.conj(matrixQ2)), data_input))**2)
        iteration_radius = [0 for i in range(tx_dimension)]
        '''
        iteration_radius[i] = np.sqrt(radius**2 - np.linalg.norm(data_input)**2 +
                np.linalg.norm(np.matmul(channel,tx_symbol_estimate))**2)

        lower_bound = np.ceil(np.real(iteration_estimate-(iteration_radius[i]/matrixR[i][i]))) + np.ceil(
                np.imag(iteration_estimate-(iteration_radius[i]/matrixR[i][i])))*1j

        upper_bound = np.floor(np.real(iteration_estimate+(iteration_radius[i]/matrixR[i][i]))) + np.floor(
                np.imag(iteration_estimate+(iteration_radius[i]/matrixR[i][i])))*1j
        '''
        iteration_radius[i] = np.sqrt(radius**2 - 
                np.linalg.norm(np.matmul(np.conj(np.transpose(matrixQ2)),data_input))**2)

        lower_bound = np.ceil(np.real((iter_y-iteration_radius[i])/matrixR[i][i])) + np.ceil(
                np.imag((iter_y-iteration_radius[i])/matrixR[i][i]))

        upper_bound = np.floor(np.real((iter_y+iteration_radius[i])/matrixR[i][i])) + np.floor(
                np.imag((iter_y+iteration_radius[i])/matrixR[i][i]))
        
        #estimate = [0 for i in range(tx_dimension)]
        #distance = [0 for i in range(tx_dimension)]


        estimate = [0 for i in range(tx_dimension)]
        estimate[i] = lower_bound - 1

        final_estimate = [0 for i in range(tx_dimension)]

        counter = 0 
        while 1: #i < tx_dimension:
            #Step 3
            print(i, counter)
            counter+=1
            estimate[i] += 1
            if estimate[i] <= upper_bound:
                #Step 5
                if i==0:
                    #Step 6
                    final_estimate[i] == iteration_estimate


                else:
                    #Remaining of step 5
                    #old_estimate = iteration_estimate
                    i -= 1

                    #y[i] = y[] - sum()

                    iteration_radius[i] = np.sqrt(iteration_radius[i+1]**2 - (y[i+2] - r[i+1][i+1]*estimate[i+1])**2)

                    iteration_estimate = tx_symbol_estimate[i] - sum(
                            (estimate[j] - tx_symbol_estimate[j])*matrixR[i][j]/matrixR[i][i]
                            for j in range(i+1,tx_dimension))

                    lower_bound = np.ceil(np.real((iter_y-iteration_radius[i])/matrixR[i][i])) + np.ceil(
                            np.imag((iter_y-iteration_radius[i])/matrixR[i][i]))

                    upper_bound = np.floor(np.real((iter_y+iteration_radius[i])/matrixR[i][i])) + np.floor(
                            np.imag((iter_y+iteration_radius[i])/matrixR[i][i]))
                    
                    estimate[i] = lower_bound - 1

                    '''
                    iteration_radius[i] = iteration_radius[i+1]-(matrixR[i+1][i+1]**2)*(estimate[i+1]-old_estimate)

                    lower_bound = np.ceil(np.real(iteration_estimate-(iteration_radius[i]/matrixR[i][i]))) + np.ceil(
                            np.imag(iteration_estimate-(iteration_radius[i]/matrixR[i][i])))*1j

                    upper_bound = np.floor(np.real(iteration_estimate+(iteration_radius[i]/matrixR[i][i]))) + np.floor(
                        np.imag(iteration_estimate+(iteration_radius[i]/matrixR[i][i])))*1j

                    estimate[i] = lower_bound - 1

                    '''

            else:
                #Step 4
                i += 1

        return estimate
