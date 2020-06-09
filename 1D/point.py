import numpy as np
import math
import matplotlib.pyplot as plt

def H(J,h,d):
    return -2*J*(d**2)-h*d

def F(J,h,T,d):
    return -2*J*(d**2) -h*d + T*((1+d)*math.log1p(d)+(1-d)*math.log1p(-d))/2-T*math.log(2)

def FXvT(J):
    x = np.linspace(-1,1,100)
    d = x[1:-1]
    Temp = [0,.5,1,2,4]
    
    vF=np.vectorize(F)
    
    for T in Temp:
        L="T="+str(T)
        Y= vF(J,1,T,d)
        Y=np.insert(Y,0,H(J,1,-1))
        Y=np.append(Y,H(J,1,1))
        plt.plot(x,Y,label=L)

    plt.ylabel('Free Energy')
    plt.xlabel('Delta')
    
    title="J="+str(J)
    plt.title(title)
    plt.legend()
    
    plt.show()

def FvHT(J):
    return 0
