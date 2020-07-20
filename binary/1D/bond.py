import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

def xlx(x):
    if x==0:return 0
    return x*math.log(x)

def H(J,h,x,y):
    return -2*J*y-h*x

def S(x, y):
    Sx,Sy=0,0
    if abs(x-.5)<0.5:
        Sx=xlx(x)+xlx(1-x)
        if y<=x and y<=1-x: 
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
    return -1 

def F(J,h,x,y,T):
    if T==0: return H(J,h,x,y)
    return H(J, h, x, y)-T*S(x,y)

def Test():
    rranges = (slice(0,1,0.01), slice(0,0.5,0.01))
    res=opt.brute(Sp, rranges, full_output=True, finish=opt.fmin)
    return res

def FvT(J, hpj=1):
    Temp = [0,.25,.5,1,2,3,6,10]
    for Tpj in Temp:
        fun = lambda x,y: F(J, hpj*abs(J), x, y, Tpj*abs(J))


def XYvT(J, hpj=1):
    Temp = np.linspace(0,10,201) 
    rranges=(slice(0,1,0.01), slice(0,0.5,0.01))
    
    con = opt.LinearConstraint([[1,0],[-1,1],[1,1]],[0,-np.inf,-np.inf],[1,0,1])
    bound = opt.Bounds([0,0],[1,0.5])
    
    mX, mY=[],[]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
       
    for Tpj in Temp:
        fun = lambda x: F(J, hpj*abs(J), x[0], x[1], Tpj*abs(J))
        res=opt.minimize(fun, [0.5,0.25], method='trust-constr', constraints=con, bounds=bound)
        mX.append(res.x[0])
        mY.append(res.x[1])
        working="Calculating T/J:"+str(Tpj)
        print(working)

    ax.plot(mX, mY, Temp, label='min')
    ax.set_xlabel('Composition: X')
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_zlabel('Temperature (T/J)')

    plt.show()

