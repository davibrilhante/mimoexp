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
                '''
                red_real = self.lllRealReduction(np.real(ch).T, delta).T
                red_imag = self.lllRealReduction(np.imag(ch).T, delta).T
                h_red.append(red_real + 1j*red_imag)
                '''

                # the real valued lll needs a lattice conversion from complex to real
                # [[re(H) -im(H)] [im(H) re(H)]]
                real = self.complex2real(ch.T)

                # The real valued lll has a different delta
                delta_real = delta - 0.25

                # The real reduction over the real converted lattice
                red_real = self.lllRealReduction(real,delta_real)

                # Now, the reduced basis is converted back to comples
                h_red.append(self.real2complex(red_real).T)

            else:
                raise NotImplementedError

            h_inv = np.linalg.pinv(ch)
            U.append(np.matmul(h_inv, h_red[-1]))

        return np.array(U), np.array(h_red)

    def complex2real(self,basis):
        rows,columns = basis.shape
        basis_real = np.zeros([2*rows,2*columns])
        for i in range(rows):
            for j in range(columns):
                basis_real[i*2:(i+1)*2,j*2:(j+1)*2] = np.array([np.real(basis[i,j]), 
                    -1*np.imag(basis[i,j]), np.imag(basis[i,j]), np.real(basis[i,j])]).reshape([2,2])
        return basis_real

    def real2complex(self,basis):
        rows, columns = basis.shape
        basis_reduced = np.zeros([rows//2,columns//2],dtype=complex)
        for i in range(rows//2):
            for j in range(columns//2):
                basis_reduced[i,j] = basis[2*i,2*j] - 1j*basis[2*i,(2*j)+1]

        return basis_reduced


    def lllComplexReduction(self, basis, delta):
        #Each column is a base vector, so we make M the number of columns
        M = np.shape(basis)[1]

        #QR orthogonalization via Gram-Schmidt
        Q, R = self.GramSchmidt(basis, norm=True)
        #Q,R = np.linalg.qr(basis)

        #Squared length of orthogonal vectors
        Rdiag = np.diag(R)
        beta = abs(Rdiag)**2

        Mu = R/(np.diag(Rdiag)@np.ones((M,M)))
        Mu = Mu.T
        #print(Mu)

        k = 1
        i_iteration = -1

        while (i_iteration < 100*M**2):
            i_iteration += 1

            #print(round(abs(np.real(Mu[k,k-1])),5))
            #print(abs(np.imag(Mu[k,k-1])))
            if (round(abs(np.real(Mu[k,k-1])),5) > 0.5 or round(abs(np.imag(Mu[k,k-1])),5) > 0.5):
                basis, Mu = self.size_reduce(basis,Mu,k,k-1)
                #print(basis)
                #print(Mu)


            #print(beta[k])
            #print((delta - abs(Mu[k,k-1])**2)*beta[k-1])
            #print()
            if (beta[k] < (delta - abs(Mu[k,k-1])**2)*beta[k-1]):
                #Swap k with k-1 if k-1th vector is longer than kth
                b = np.copy(basis[:,k])
                basis[:,k] = basis[:,k-1]
                basis[:,k-1] = b
                #print(basis[:,k])
                #print(basis[:,k-1])

                #Swap the kth column and the k-1th column
                muswap = np.copy(Mu[k-1,0:k-1])
                Mu[k-1,0:k-1] = Mu[k,0:k-1]
                Mu[k,0:k-1] = muswap

                #The might be in these indexes
                old_muk = np.copy(Mu[k+1:M,k])
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
                #The error might be here
                #print(basis)
                for i in range(k-1,0,-1):
                    #print(abs(np.real(Mu[k,k-1])))
                    #print(abs(np.imag(Mu[k,k-1])))
                    if (round(abs(np.real(Mu[k,k-1])),5) > 0.5 or round(abs(np.imag(Mu[k,k-1])),5) > 0.5):
                        #print('Entrou')
                        basis, Mu = self.size_reduce(basis,Mu,k,i)

                if k<M-1:
                    k += 1

                else:
                    basis_reduced = np.copy(basis)
                    Mu = abs(Mu)
                    return basis_reduced

        # i_iteration exceeded the limit, so return basis so far.
        # This is the sub-optimal solution
        basis_reduced = np.copy(basis)
        return basis_reduced


    def lllRealReduction(self, basis, delta):
        basis = np.array(basis)
        Q, R = self.GramSchmidt(basis)

        m, n = np.shape(basis)
        k=1
        counter = 0
        while k < n:
            #print('===========',k,'===========')
            for j in range(k-1, -1, -1):
                mu = self.mu(basis, Q, k, j)
                #print('mu_{k}_{j}={mu}'.format(k=k,j=j,mu=mu))

                #check whether reduction is required
                if (abs(np.real(mu))>0.5 or abs(np.imag(mu))>0.5):
                    basis[:,k] = basis[:,k] - basis[:,j]*round(mu)
                    Q, R = self.GramSchmidt(basis)

            # Recalculate mu for the lovasz condition below
            mu = self.mu(basis, Q, k, k-1)
            #Lovasz condition
            if Q[:,k].conj().T@Q[:,k] >= (delta - mu**2)*(Q[:,k-1].conj().T@Q[:,k-1]):
                #print('Lovasz Ok')
                k = k+1
            else:
                #print('Needs Swap')
                old_b = np.copy(basis[:,k])
                basis[:,k] = np.copy(basis[:,k-1])
                basis[:,k-1] = np.copy(old_b)

                Q, R = self.GramSchmidt(basis)
                k = max(k-1,1)

        return basis

    def mu(self, B, Q, i, j):
        v = B[:,i]
        u = Q[:,j]
        #print(u,v)
        return (u.conj().T@v)/(u.conj().T@u)

    def size_reduce(self, basis, Mu, k, j):
        eta = round(Mu[k,j])
        basis[:,k] = basis[:,k] - eta*basis[:,j]
        #print(basis)

        # Perhaps the error is here!
        for i in range(j):
            #print(Mu[k,i], Mu[j,i])
            Mu[k,i] = Mu[k,i] - eta*Mu[j,i]

        #print(Mu)

        Mu[k,j] = Mu[k,j] - eta

        return basis, Mu


    def GramSchmidt(self, basis, norm=False):
        m, n = np.shape(basis)

        # Q is a  MxN Orthogonal matrix and R is a NxN upper triangular matrix
        Q = np.zeros([m,n]) 
        R = np.zeros([n,n])


        if type(basis[0,0]) == np.complex128:
            Q = np.zeros([m,n], dtype=type(basis[0,0]))
            R = np.zeros([n,n], dtype=type(basis[0,0]))


        # Q[n] is equal the vector basis[n] minus the projection on the other vectors
        for i in range(n):
            v = basis[:,i]
            
            #for j in range(i - 1):
            for j in range(i):
                # R is the value of the projection of vector Q[i] on vector Q[j]
                #R[i][j] = (Q[j].conj().T@Q[i])/np.inner(Q[j],Q[j])
                R[j,i] = (basis[:,i]@Q[:,j].conj().T)/np.inner(Q[:,j],Q[:,j])

                # Orthogonalization
                v = v - R[j,i]*Q[:,j]

            # Normalization is optional
            if norm:
                v_norm = np.linalg.norm(v)
                Q[:,i] = v/v_norm
                R[i,i] = v_norm

            else:
                Q[:,i] = v
                R[i,i] = 1

        return Q, R
