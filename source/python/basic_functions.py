import numpy as np 
import math
import matplotlib.pyplot as plt

from awgn import *

#basic functions
def gen_bin(K, ns=128):
  #K = number of terminals
  #ns = stream size
  return np.random.randint(2, size=(K,ns))

def mod_bpsk(x):
  qq=x +0j
  qq[x==0]=-1 
  qq[x==1]=1 
  return qq

def demod_bpsk(q_received):
   xx=q_received
   xx[xx.real<0]=0+0j
   xx[xx.real>0]=1+0j
   return xx

def ber(x_tx, x_rx):
  l = len(x_tx)
  errors = sum(x_tx!=x_rx)
  return errors/l

def ber_average_per_user(X_tx, X_rx):
    KK = X_tx.shape[0]
    bersum = 0
    for i in range(KK):
      bersum = ber(X_tx[i,:],X_rx[i,:])+bersum
    return bersum/KK

def rayleigh_channel(M,K):
    #TRUE CHANNEL 
    Gg =  np.random.normal(0,1,(M,K))+1j*(np.random.normal(0,1,(M,K)))        
    for k in range(K):
        Gg[:,k] = math.sqrt(1/np.var(Gg[:,k]))*Gg[:,k]
    return Gg


def transmitted_message(q_pred,G,snr):
  #return (np.dot(np.transpose(G),q_pred)+np.random.normal(0,10**(-snr),np.dot(np.transpose(G),q_pred).shape))
   return awgn(np.dot(np.transpose(G),q_pred),snr)
def zf(G):
  #Zero-forcing preconding
  return  np.dot(np.conj(G),np.linalg.inv(np.dot(np.transpose(G),np.conj(G))));

def sending_stream(M,K,Ns,G,P,snr):
  x = gen_bin(K,Ns)
  q = mod_bpsk(x)
  q_pred = np.dot(P,q)
  y = transmitted_message(q_pred,G,snr)
  x_demod = demod_bpsk(y)
  return (ber_average_per_user(x,x_demod))


#  def channel_estimation(G,SNR,Np,mode):
#      K,M = G.shape
#      pilots = gen_bin(K,Np)  
#      beta_k = (1./(10**(SNR/10)))
