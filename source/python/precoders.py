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
            red_real = self.lllRealReduction(np.real(ch).T, delta).T
            red_imag = self.lllRealReduction(np.imag(ch).T, delta).T

            h_red.append(red_real + 1j*red_imag)
            h_inv = np.linalg.pinv(ch)
            U.append(np.matmul(h_inv, h_red[-1]))



        return np.array(U), np.array(h_red)


    def lllComplexReduction(self, basis, delta):
        M = np.shape(basis)[1]
        Q, R = self.GramSchmidt(basis.T)
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
        Q, R = self.GramSchmidt(basis)

        m, n = np.shape(basis)
        k=1
        counter = 0
        while k < n:
            for j in range(k-1, -1, -1):
                mu = self.mu(basis, Q, k, j)
                #print(k, j, mu)

                #check whether reduction is required
                if (abs(np.real(mu))>0.5 or abs(np.imag(mu))>0.5):
                    basis[:,k] = basis[:,k] - basis[:,j]*round(mu)
                    Q, R = self.GramSchmidt(basis)

            mu = self.mu(basis, Q, k, k-1)
            #Lovasz condition
            if Q[:,k].conj().T@Q[:,k] >= (delta - mu**2)*(Q[:,k-1].conj().T@Q[:,k-1]):
                k = k+1
            else:
                temp = np.copy(basis[:,k])
                basis[:,k] = np.copy(basis[:,k-1])
                basis[:,k-1] = np.copy(temp)

                Q, R = self.GramSchmidt(basis)
                k = max(k-1,1)

        return basis

    def mu(self, B, Q, i, j):
        v = B[:,i]
        u = Q[:,j]
        return (u.conj().T@v)/(u.conj().T@u)

    def GramSchmidt(self, basis, norm=False):
        m, n = np.shape(basis)

        # Q is a  MxN Orthogonal matrix and R is a NxN upper triangular matrix
        Q = np.zeros([m,n])
        R = np.zeros([n,n])

        # Q[n] is equal the vector basis[n] minus the projection on the other vectors
        for i in range(n):
            v = basis[:,i]
            
            #for j in range(i - 1):
            for j in range(i):
                # R is the value of the projection of vector Q[i] on vector Q[j]
                #R[i][j] = (Q[j].conj().T@Q[i])/np.inner(Q[j],Q[j])
                R[j][i] = (basis[:,i]@Q[:,j].conj().T)/np.inner(Q[:,j],Q[:,j])

                # Orthogonalization
                v = v - R[j][i]*Q[:,j]

                #Q[:,i] = Q[:,i] - R[i][j]*Q[:,j]
                # Q[i] is now orthogonal
                
            # Normalization is optional
            if norm:
                v_norm = np.linalg.norm(v)
                Q[:,i] = v/v_norm
                R[i][i] = v_norm

            else:
                Q[:,i] = v
                R[i][i] = 1

        return Q, R
