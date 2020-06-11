import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def H(J,h,x):
    return -J*x*(1-x)-h*x

def S(J,h,x):
    return -x*math.log(x)-(1-x)*math.log(1-x)

def F(J,h,T,x):
    if abs(x-.5)==.5:
        return H(J,h,x)
    return H(J,h,x)-T*S(J,h,x)

def FvXT(J, hpj=1):
    x = np.linspace(0,1,1000)
    Temp = np.linspace(0, 10, 51) 
    
    Mx=[]
    MT=[]
    

    vF=np.vectorize(F)
    
    show=0
    
    for Tpj in Temp:
        L="T/J="+str(Tpj)
        Y= vF(J,hpj*abs(J),Tpj*abs(J),x)
        if show%5==0:
            plt.plot(x,Y,label=L)

        fun = lambda x: F(J,h,T,x)
        x0=0.5
        xmin=minimize(fun, x0, method='Nelder-Mead', tol=1e-9)
        
        show+=1

    plt.ylabel('Free Energy')
    plt.xlabel('Delta')
    
    title="J="+str(J)+", h/J="+str(hpj)
    plt.title(title)
    plt.legend(loc='center left', bbox_to_anchor=(0.8,0.5))
    
    plt.show()

def Min():
    J,h,T=1,1,1
    return xmin.x
