import numpy as np
import math
import matplotlib.pyplot as plt

def H(J,h,d):
    return -2*J*(d**2)-h*d

def S(J,h,d):
    return -((1+d)*math.log1p(d)+(1-d)*math.log1p(-d))/2+math.log(2)

def F(J,h,T,d):
    if abs(d)==1:
        return H(J,h,d)
    else:
        return H(J,h,d)-T*S(J,h,d)

def FXvT(J, hpj=1):
    d = np.linspace(-1,1,100)
    Temp = [0,.5,1,2,3,4,6,10]
    
    vF=np.vectorize(F)
    
    for Tpj in Temp:
        L="J/T="+str(Tpj)
        Y= vF(J,hpj*abs(J),Tpj*abs(J),d)
        plt.plot(d,Y,label=L)

    plt.ylabel('Free Energy')
    plt.xlabel('Delta')
    
    title="J="+str(J)+", h/J="+str(hpj)
    plt.title(title)
    plt.legend()
    
    plt.show()

def FvHT(J):
    return 0
