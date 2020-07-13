import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

Z=4

def xlx(x):
    if x<0:return -1
    if x==0:return 0
    if x>0.5:
        return x*math.log1p(x-1)
    return x*math.log(x)

def H(J,h,x,y):
    #return -Z*J*y-h*x
    return -J*(Z/2)*(4*y-1)-h*(2*x-1)

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

def min(J=1, hpj=1, samp=200, Trang=[0,4]):
    Temp = np.linspace(Trang[0],Trang[1],samp+1) 
    
    con = opt.LinearConstraint([[-1,1],[1,1]],[-np.inf,-np.inf],[0.0,1.0])
    bound = opt.Bounds([0.0,0.0],[1.0,0.5])
    guess = [0.,0.]
    
    Free, E, C, mX, mY=[],[],[],[],[]
    delta = (Trang[1]-Trang[0])/samp
    for i in range(len(Temp)):
        working="Calculating: T/J="+str(Temp[i])
        print(working)
        fun = lambda x: F(J, hpj*abs(J), x[0], x[1], Temp[i]*abs(J))
        res=opt.minimize(fun, guess, method='trust-constr', constraints=con, bounds=bound)
        mX.append(res.x[0])
        mY.append(res.x[1])
        Free.append(fun(res.x))
        guess=res.x
        if i>1:
            #Calculating E
            dF=Free[i]-Free[i-2]
            dF/=(2*delta)
            E.append(Free[i-1]-Temp[i-1]*dF)
            #Calculating C
            ddF=Free[i]-2*Free[i-1]+Free[i-2]
            ddF/=(delta**2)
            C.append(-Temp[i-1]*ddF)

    fig = plt.figure(figsize=plt.figaspect(2.))
    fig.suptitle('2D square lattice, bond approx')
    
    #plot min free energy wrt temp
    print('Plotting free energy')
    ax = fig.add_subplot(2,2,1)
    ax.set_xlabel('Temperature (T/J)')
    ax.set_ylabel('Min Free Energy')
    ax.plot(Temp, Free, label='min free energy')
   
 
    #plot composition wrt temp
    print('Plotting composition')
    ax=fig.add_subplot(2,2,2, projection='3d')
    ax.plot(mX, mY, Temp, label='min')
    ax.set_xlabel('Composition: X')
    ax.set_xlim(0.0,1.0)
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_ylim(0.0,0.5)
    ax.set_zlabel('Temperature (T/J)')
    ax.set_zlim(0,Trang[1])

    #plot E wrt temp
    Tp=np.linspace(Trang[0]+delta, Trang[1]-delta, samp-1)
    ax=fig.add_subplot(2,2,3)
    ax.plot(Tp, E)
    ax.set_xlabel('Temperature (T/J)')
    ax.set_ylim(-5,5)
    ax.set_ylabel('E=F-T*dF/dT')

    #plot C wrt temp
    ax=fig.add_subplot(2,2,4)
    ax.plot(Tp, C)
    ax.set_xlabel('Temperature (T/J)')
    ax.set_ylabel('C=-T*d2F/dT2')
    ax.set_ylim(-5,5)
    
    #high temp slope, should be ln(2)~0.693
    slope = Free[samp]-Free[math.floor(samp*0.75)]
    slope /= (Temp[samp]-Temp[math.floor(samp*0.75)])
    print('High temp slope: ' + str(slope))
    #intercept
    intercept=Free[samp]-slope*Temp[samp]
    print('Intercept: ' + str(intercept))

    #display plots
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
        #while j<=1: 
        while j<=i and j<=(1-i):
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
    ax.set_xlim(-0.1,1.1)
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_ylim(-0.1,0.6)
    ax.set_zlabel('Free Energy')
    plt.show()
