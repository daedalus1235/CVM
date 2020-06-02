import math
import numpy as np
import matplotlib.pyplot as plt
import time

def Sig(N):
    sig=[]
    for x in range(2**N):
        sigi=[]
        i=0
        while i<N:
            if x%2==0:
                sigi.append(-1)
            else:
                sigi.append(1)
            x=math.floor(x/2)
            i+=1
        sig.append(sigi)
    return sig

def Part(N, B, h, J):
    Z=0
    for s in Sig(N):
        Z+=math.exp(-B*Ham(s,h,J))
    return Z

def Ham(s, h, J):
    Hh=0
    Hj=0
    for spin in s:
        Hh+=spin
    for x in range(len(s)-1):
        Hj+=s[x]*s[x+1]
    H=h*Hh+J*Hj
    return H

def Ent(s, h, J):
    return 0


def minH(N,h,J):
    start = time.time()
    configs=Sig(N)
    Hmin=Ham(configs[0],h,J)
    imin=[]
    i=0

    for state in configs:
        H=Ham(state, h, J)
        if H==Hmin:
            imin.append(i)
        if H<Hmin:
            Hmin=H
            imin=[]
            imin.append(i)
        i+=1
    
    print("Energy: ", Hmin)
    n=0
    for state in imin:
        n+=1
        print("State ", n,": ", configs[state])
    
    print("Execution time:", time.time()-start, "s")
    
    return Hmin
