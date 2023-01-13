#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         
import numpy as np
from detectors import ZeroForcingSic
import sage.all
from sage.all import Matrix,vector,ZZ, QQ, CC
from olll import reduction

class Precoder:
    def __init__(self):
        pass

    def precode(self, *args, **kwargs):
        raise NotImplementedError

class LatticeReductionAided(Precoder):
    def __init__(self):
        super().__init__()

    def precode(self, channel, delta=3/4):
        h_red = []
        U = []

        if not isinstance(channel,list):
            channel = [channel]

        for ch in channel:
            #h_1 = Matrix(QQ,np.matrix(np.real(ch))).LLL(delta)
            #h_1 = reduction(np.real(ch), delta)
            #print(h_1)
            red_real = self.lllRealReduction(np.real(ch), delta)
            red_imag = self.lllRealReduction(np.imag(ch), delta)

            h_red.append(red_real + 1j*red_imag)
            h_inv = np.linalg.pinv(ch)
            U.append(np.matmul(h_inv, h_red[-1]))



        return np.array(U), np.array(h_red)


    def lllComplexReduction(self, basis, delta):
        M = np.shape(basis)[1]
        Q, R = self.GramSchmidt(basis)
        print(basis)
        print(Q@R)

        '''
        beta = abs()
        mu = 

        k=1
        i_counter = 0
        '''



    def lllRealReduction(self, basis, delta):
        basis = np.array(basis)
        Q = self.GramSchmidt(basis)

        m,n = basis.shape
        k=1
        counter = 0
        while k < n:
            for j in range(k-1, -1, -1):
                mu = self.mu(basis, Q, k, j)

                #check whether reduction is required
                if (abs(np.real(mu))>0.5 or abs(np.imag(mu))>0.5):
                    basis[k] = basis[k] - basis[j]*round(mu)
                    Q = self.GramSchmidt(basis)

            mu = self.mu(basis, Q, k, k-1)
            #Lovasz condition
            if Q[k].conj().T@Q[k] >= (delta - mu**2)*(Q[k-1].conj().T@Q[k-1]):
                k = k+1
            else:
                temp = np.copy(basis[k])
                basis[k] = np.copy(basis[k-1])
                basis[k-1] = np.copy(temp)

                Q = self.GramSchmidt(basis)
                k = max(k-1,1)

        return basis

    def mu(self, B, Q, i, j):
        v = B[i]
        u = Q[j]
        return (u.conj().T@v)/(u.conj().T@u)

    def GramSchmidt(self, basis):
        m, n = np.shape(basis)

        # Q is a  MxN Orthogonal matrix and R is a NxN upper triangular matrix
        Q = np.zeros([m,n])
        R = np.zeros([n,n])

        # Q[n] is equal the vector basis[n] minus the projection on the other vectors
        #Q[0] = basis[0]
        #for i in range(1,n):
        for i in range(n):
            Q[i] = basis[i]
            
            for j in range(i):
                # R is the value of the projection if vector Q[i] on vector Q[j]
                #R[i][j] = (Q[j].conj().T@Q[i])/np.inner(Q[j],Q[j])
                R[i][j] = (basis[i]@Q[j].conj().T)/np.inner(Q[j],Q[j])

                # Q[i] is now orthogonal
                Q[i] = Q[i] - R[i][j]*Q[j]
                
            # Now Q[i] is orthonormal
            Q[i] = Q[i]/np.linalg.norm(Q[i])

            # R is equal to the dot product between normalized Q[i] and basis[i]
            R[i][i] = np.matrix(Q[i])@basis[i].conj().T#/np.linalg.norm(Q[i])/np.linalg.norm(Q[i])


        return Q, R.conj().transpose()
