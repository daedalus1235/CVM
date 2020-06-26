import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

Z=4

def xlx(x):
    if x==0:return 0
    return x*math.log(x)

def H(J,h,x,y):
    return -Z*J*y-h*x

def S(x, y):
    Sx,Sy=0,0
    if abs(x-.5)<0.5:
        Sx=xlx(x)+xlx(1-x)
        Sx*=(Z-1)
        if y<=x and y<=1-x: 
            Sy=xlx(x-y)+2*xlx(y)+xlx(1-x-y)
            Sy*=(Z/2)
            return Sx-Sy
    return -1 

def F(J,h,x,y,T):
    if T==0: return H(J,h,x,y)
    return H(J, h, x, y)-T*S(x,y)

def XYvT(J, hpj=1):
    Temp = np.linspace(0,10,201) 
    rranges=(slice(0,1,0.01), slice(0,0.5,0.01))
    
    con = opt.LinearConstraint([[1,0],[-1,1],[1,1]],[0,-np.inf,-np.inf],[1,0,1])
    bound = opt.Bounds([0,0],[1,0.5])
    
    Free, mX, mY=[],[],[]
    
    for Tpj in Temp:
        fun = lambda x: F(J, hpj*abs(J), x[0], x[1], Tpj*abs(J))
        res=opt.minimize(fun, [0.5,0.25], method='trust-constr', constraints=con, bounds=bound)
        mX.append(res.x[0])
        mY.append(res.x[1])
        Free.append(fun(res.x))
        working="Calculating T/J:"+str(Tpj)
        print(working)
    
    fig = plt.figure(figsize=plt.figaspect(2.))
    fig.suptitle('2D square lattice, bond approx')
    
    #plot min free energy wrt temp
    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('Temperature (T/J)')
    ax.set_ylabel('Min Free Energy')
    ax.plot(Temp, Free, label='min free energy')
   
    
    #plot composition wrt temp
    ax=fig.add_subplot(2,1,2, projection='3d')
    ax.plot(mX, mY, Temp, label='min')
    ax.set_xlabel('Composition: X')
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_zlabel('Temperature (T/J)')

    plt.show()

