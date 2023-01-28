#! /usr/bin/env python3                                                         
# -*- docing : utf8 -*-                                                         
import numpy as np
from detectors import ZeroForcingSic
#import sage.all
#from sage.all import Matrix,vector,ZZ, QQ, CC
#from olll import reduction

class Precoder:
    def __init__(self):
        pass

    def precode(self, *args, **kwargs):
        raise NotImplementedError

class LatticeReductionAided(Precoder):
    def __init__(self):
        super().__init__()
        self.type = 'real'

    def precode(self, channel, delta=3/4):
        h_red = []
        U = []

        if not isinstance(channel,list):
            channel = [channel]

        for ch in channel:
            #h_1 = Matrix(QQ,np.matrix(np.real(ch))).LLL(delta)
            #h_1 = reduction(np.real(ch), delta)
            #print(h_1)
            if self.type == 'complex':
                h_red.append(self.lllComplexReduction(ch,delta))

            elif self.type == 'real':
                red_real = self.lllRealReduction(np.real(ch).T, delta).T
                red_imag = self.lllRealReduction(np.imag(ch).T, delta).T
                h_red.append(red_real + 1j*red_imag)

            else:
                raise NotImplementedError

            h_inv = np.linalg.pinv(ch)
            U.append(np.matmul(h_inv, h_red[-1]))

        return np.array(U), np.array(h_red)


    def lllComplexReduction(self, basis, delta):
        #Each column is a base vector, so we make M the number of columns
        M = np.shape(basis)[1]

        #QR orthogonalization via Gram-Schmidt
        #Q, R = self.GramSchmidt(basis)
        Q,R = np.linalg.qr(basis)

        #Squared length of orthogonal vectors
        beta = abs(np.diag(R))**2

        Mu = R/(np.diag(np.diag(R))@np.ones(M))
        Mu = Mu.T

        k = 1
        i_iteration = 0

        while (i_iteration < 100*M**2):
            i_iteration += 1
            if (abs(np.real(Mu[k,k-1])) > 0.5) or (abs(np.imag(Mu[k,k-1])) > 0.5):
                basis, Mu = self.size_reduce(basis,Mu,k,k-1)


            if (beta[k] < delta - beta[k-1]*abs(Mu[k,k-1])**2):
                #Swap k with k-1 if k-1th vector is longer than kth
                b = basis[:,k]
                basis[:,k] = basis[:,k-1]
                basis[:,k-1] = b

                #Swap the kth column and the k-1th column
                muswap = Mu[k-1,:k-2]
                Mu[k-1,:k-2] = Mu[k,:k-2]
                Mu[k,:k-2] = muswap

                old_muk = Mu[k+1:M,k]
                old_beta1 = beta[k-1]
                old_betak = beta[k]
                old_mu = Mu[k,k-1]

                Mu[k+1:M,k] = Mu[k+1:M,k-1] - Mu[k+1:M,k]*Mu[k,k-1]
                beta[k-1] = beta[k] + beta[k-1]*abs(Mu[k,k-1])**2
                beta[k] = beta[k]*old_beta1/beta[k-1]
                Mu[k,k-1] = Mu[k,k-1].T*old_beta1/beta[k-1]
                Mu[k+1:M,k-1] = Mu[k+1:M,k-1]*Mu[k,k-1] + old_muk*old_betak/beta[k-1]

                if k>1:
                    k -= 1

            else:
                for i in range(k-3,0,-1):
                    if (abs(np.real(Mu[k,k-1])) > 0.5 or abs(np.imag(Mu[k,k-1])) > 0.5):
                        basis, Mu = self.size_reduce(basis,Mu,k,i)
                if k<M-1:
                    k += 1
                else:
                    basis_reduced = basis
                    Mu = abs(Mu)
                    return basis_reduced

        basis_reduced = basis
        return basis_reduced


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

    def size_reduce(self, basis, Mu, k, j):
        eta = round(Mu[k,j])
        basis[:,k] = basis[:,k] - eta*basis[:,j]
        for i in range(j-1):
            Mu[k,i] = Mu[k,i] - eta*Mu[j,i]

        Mu[k,j] = Mu[k,j] - eta

        return basis, Mu


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
