import numpy as np
import math
import matplotlib.pyplot as plt

def H(J,h,d):
    return -2*J*(d**2)-h*d

def F(J,h,T,d):
    return -2*J*(d**2) -h*d + T*((1+d)*math.log1p(d)+(1-d)*math.log1p(-d))/2+T*math.log(2)

def Plot(J):
    d = np.linspace(-1,1,100)
    d = d[1:-1]
    Temp = [0,.5,1,2,4]
    
    vF=np.vectorize(F)
    
    for T in Temp:
        L="T="+str(T)
        plt.plot(d, vF(J,1,T,d), label=L)

    plt.ylabel('Free Energy')
    plt.xlabel('Delta')
    
    title="J="+str(J)
    plt.title(title)
    plt.legend()
    
    plt.show()
