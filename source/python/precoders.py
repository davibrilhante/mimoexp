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
        basis = np.array(basis)

        m, n = basis.shape

        norm_basis = np.zeros([n])

        for j in range(n):
            norm_basis[j] = np.inner(np.conjugate(basis[:,j]),basis[:,j])

        mu = np.zeros([n,n])

        for j in range(n):
            for i in range(j+1,n):
                inner = np.inner(np.conjugate(basis[:,j]),basis[:,i])
                facsum = sum([np.conjugate(mu[j][k])*mu[i][k]*norm_basis[k] for k in range(1,j-1)])
                mu[i][j] = (1/norm_basis[j])*(inner - facsum)
                norm_basis[i] = norm_basis[i] - norm_basis[j]
            #endfor
        #endfor

        U = np.identity(n)
        k = 1

        while (k <= n-1):
            if (np.real(mu[k][k-1])>0.5) or (np.imag(mu[k][k-1])>0.5):
                basis, U, mu = self.sizeReduce(basis,U,mu,k,k-1)

            if (norm_basis[k] >= 
                (delta - (mu[k][k-1]*np.conjugate(mu[k][k-1])))*norm_basis[k-1]):
                basis, norm_basis, mu = self.swapUpdate(basis, norm_basis, mu, k)
                
                temp = U[:,k]
                U[:,k] = U[:,k-1]
                U[:,k-1] = temp

                k = max(1,k-1)
            else:
                for j in reversed(range(k-1)):
                    if (np.real(mu[k][j])>0.5) or (np.imag(mu[k][j])>0.5):
                        basis, U, mu = self.sizeReduce(basis,U,mu,k,j)
                k = k+1

        return basis, U


    def sizeReduce(self,basis, U, mu, k, j):
        c = int(np.real(mu[k][j])) + 1j*int(np.imag(mu[k][j]))
        basis[:,k] = basis[:,k] - c*basis[:,j]
        U[:,k] = U[:,k] - c*U[:,j]

        for l in range(j):
            mu[k][l] = mu[k][l] - c*mu[j][l]

        return basis, U, mu



    def swapUpdate(self,basis,norm_basis, mu, k):
        m, n = basis.shape
        new_basis_k_1 = basis[:,k-1]
        basis[:,k-1] = basis[:,k]
        basis[:,k] = new_basis_k_1

        norm_copy = np.copy(norm_basis)
        mu_copy = np.copy(mu)

        norm_mu = mu[k][k-1]*np.conjugate(mu[k][k-1])
        norm_copy[k-1] = norm_basis[k] + norm_mu*norm_basis[k-1]

        mu_copy[k][k-1] = np.conjugate(mu[k][k-1])*(norm_basis[k-1]/norm_copy[k-1])
        
        new_norm_mu = mu_copy[k][k-1]*np.conjugate(mu_copy[k][k-1])
        norm_copy[k] = norm_basis[k-1] - new_norm_mu*norm_copy[k-1]

        for i in range(k,n):
            mu_copy[i][k-1] = mu[i][k-1]*mu_copy[k][k-1] + mu[i][k]*(norm_basis[k]/norm_copy[k-1])

        for i in range(k,n):
            mu_copy[i][k] = mu[i][k-1] - mu[i][k]*mu[k][k-1]

        for j in range(0,k-1):
            mu_copy[k-1][j] = mu[k][j]

        for j in range(0,k-1):
            mu_copy[k][j] = mu[k-1][j]

        return basis, norm_copy, mu_copy


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
        Q = np.zeros([m,n])

        R = np.zeros([n,n])

        Q[0] = basis[0]
        for i in range(1,n):
            Q[i] = basis[i]
            
            for j in range(i):
                R[i][j] = (Q[j].conj().T@Q[i])/np.inner(Q[j],Q[j])
                Q[i] = Q[i] - R[i][j]*Q[j]

            #R[i][i] = Q[i].conj().T@basis[i]/np.linalg.norm(Q[i])/np.linalg.norm(Q[i])

        return Q #, R.conj().transpose()
