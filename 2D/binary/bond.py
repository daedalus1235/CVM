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
    if abs(x-.5)<=0.5:
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

def XYvT(J=1, hpj=1):
    Temp = np.linspace(0,10,201) 
    
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

def debugH(J=1,h=1):
    x,y=[],[]
    dx,dy=0.05, 0.05
    
    i=0

    while i<=1:
        j=0
        while j<=i and j<=(1-i):
            x.append(i)
            y.append(j)
            j+=dy
        i+=dx

    Ham=np.vectorize(H)

    z=Ham(J,h,x,y)

    fig = plt.figure()
    ax=fig.gca(projection='3d')

    ax.plot_trisurf(x,y,z)
    ax.set_xlabel('Composition: X')
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_zlabel('Hamiltonian')
    plt.show()

def debugTS(Temp=1):
    x,y=[],[]
    dx,dy=0.01, 0.01
    
    i=0

    while i<=1:
        j=0
        while j<=i and j<=(1-i):
            x.append(i)
            y.append(j)
            j+=dy
        i+=dx

    Ent=np.vectorize(S)

    z=Temp*Ent(x,y)

    fig = plt.figure()
    ax=fig.gca(projection='3d')

    ax.plot_trisurf(x,y,z, linewidth=0)
    ax.set_xlabel('Composition: X')
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_zlabel('Entropy')
    plt.show()


def debugF(J=1,h=1, T=1):
    x,y=[],[]
    dx,dy=0.05, 0.05
    
    i=0

    while i<=1:
        j=0
        while j<=1: #j<=i and j<=(1-i):
            x.append(i)
            y.append(j)
            j+=dy
        i+=dx

    Free=np.vectorize(F)

    z=Free(J,h,x,y,T)

    fig = plt.figure()
    ax=fig.gca(projection='3d')

    ax.plot_trisurf(x,y,z)
    ax.set_xlabel('Composition: X')
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_zlabel('Hamiltonian')
    plt.show()
