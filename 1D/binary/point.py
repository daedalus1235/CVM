import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as opt

def H(J,h,x):
    return -J*x*(1-x)-h*x

def S(x):
    if abs(x-.5)<.5:
        return -x*math.log(x)-(1-x)*math.log(1-x)
    return 0

def F(J,h,T,x):
    return H(J,h,x)-T*S(x)

def FvT(J, hpj=1):
    x = np.linspace(0,1,1000)
    Temp = np.linspace(0, 10, 51) 
    
    Mx=[]
    MF=[]
    

    vF=np.vectorize(F)
    
    show=0
    
    for Tpj in Temp:
        if show%10==0:
            Y= vF(J,hpj*abs(J),Tpj*abs(J),x)
            L="T/J="+str(Tpj)
            plt.plot(x,Y,label=L)
        
        fun = lambda x: F(J,hpj*abs(J),Tpj*abs(J),x)
        xmin=opt.fminbound(fun,0,1, xtol=1e-9)
        
        x0=xmin
        
        Mx.append(x0)
        MF.append(fun(x0))
        
        show+=1
    
    plt.plot(Mx, MF,Label='min', color='k')

    plt.ylabel('Free Energy')
    plt.xlabel('Composition')
    
    title="J="+str(J)+", h/J="+str(hpj)
    plt.title(title)
    plt.legend(loc='center left', bbox_to_anchor=(0.9,0.5))
    
    plt.show()

def XvT(J,hpj=1):
    Temp=np.linspace(0,10,51)

    xmin=[]

    for Tpj in Temp:
        fun = lambda x: F(J, hpj*abs(J), Tpj*abs(J),x)
        xmin.append(opt.fminbound(fun, 0, 1, xtol=1e-9))

    plt.plot(Temp, xmin)
    plt.ylabel("Composition")
    plt.xlabel("Temperature (T/J)")

    title="J="+str(J)+", h/J="+str(hpj)
    plt.title(title)

    plt.show()
