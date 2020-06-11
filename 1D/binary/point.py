import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as opt

def H(J,h,x):
    return -J*x*(1-x)-h*x

def S(J,h,x):
    if abs(x-.5)<.5:
        return -x*math.log(x)-(1-x)*math.log(1-x)
    return 0

def F(J,h,T,x):
    return H(J,h,x)-T*S(J,h,x)

def FvXT(J, hpj=1):
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
        x0=0.5
        xmin=opt.fminbound(fun,0,1, xtol=1e-9)
        
        x0=xmin
        
        Mx.append(x0)
        MF.append(fun(x0))
        
        show+=1
    
    plt.plot(Mx, MF,Label='min Free', color='k')

    plt.ylabel('Free Energy')
    plt.xlabel('Delta')
    
    title="J="+str(J)+", h/J="+str(hpj)
    plt.title(title)
    plt.legend(loc='center left', bbox_to_anchor=(0.9,0.5))
    
    plt.show()
