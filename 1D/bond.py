import numpy as np
import matplotlib.pyplot as plt
import math

def H(J,h,dx,y):
    return -J*(1-2*y)-h*dx

def S(dx, y):
    if abs(dx)==1:
        return 0
    #idk how to use if statements, but dy will have to be bounded by x
    Sx=((1+dx)*math.log1p(dx)+(1-dx)*math.log1p(dx))/2-math.log(2)
    Sy=((1+dx-y)*math.log1p(dx-y)+2*y*math.log(y)+(1-dx-y)*math.log1p(-dx-y))/2-math.log(2)
    return Sx-Sy

def F(J,h,dx,y,T):
    return H(J,h,dx,y)-T*S(dx,y)

def Plot(J,y=.25,hpj=1):
    dx = np.linspace(-1,1,1000)
    Temp = [0,.5,1,2,3,4,6,10]

    vF=np.vectorize(F)

    for Tpj in Temp:
        L="T/J="+str(Tpj)
        Y=vF(J,hpj*abs(J),dx,y,Tpj*abs(J))
        plt.plot(d,Y,label=L)

    plt.ylabel('Free Energy')
    plt.xlabel('Delta X')
    
    title="J="+str(J)+", h/J="+str(hpj)
    plt.title(title)
    plt.legend(loc='center', bbox_to_anchor=(1,0.5))

    plt.show()

