import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.optimize as opt

def xlx(x):
    if x==0:return 0
    return x*math.log(x)

def H(J,h,x,y):
    return -J*y-h*x

def S(x, y):
    Sx,Sy=0,0
    if abs(x-.5)<0.5:
        Sx=xlx(x)+xlx(1-x)
        if y==0:
            return Sx
        elif y<x and y<1-x: 
            Sy=xlx(x-y)+2*xlx(y)+xlx(1-x-y)
            return Sx-Sy
    return -1 

def Sp(x):
    Sx,Sy=0,0
    if abs(x[0]-.5)<=0.5:
        Sx=xlx(x[0])+xlx(1-x[0])
        if x[1]>=0 and x[1]<=x[0] and x[1]<=(1-x[0]): 
            Sy=xlx(x[0]-x[1])+2*xlx(x[1])+xlx(1-x[0]-x[1])
            return Sx-Sy
    return 1000 

def F(J,h,x,y,T):
    return H()-T*S()

def Test():
    rranges = (slice(0,1,0.01), slice(0,0.5,0.01))
    res=opt.brute(Sp, rranges, full_output=True, finish=opt.fmin)
    return res

def FvT(J, hpj=1):
    Temp = [0,.25,.5,1,2,3,6,10]
    for Tpj in Temp:
        fun = lambda x,y: F(J, hpj*abs(J), x, y, Tpj*abs(J))


def XYvT(J, hpj=1):
    Temp = [0,.25,.5,1,2,3,4,6,10]
