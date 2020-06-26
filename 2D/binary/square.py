import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

#coefficients for multiplicity
a=[1,1]         #x
b=[1,2,1]       #y
c=[1,4,4,2,4,1] #z

def xlx(x):
    if x==0: return 0
    return x*math.log(x)

def H(J,h,z):
    x = z[1]+2*z[2]+z[3]+3*z[4]+z[5]

    y = z[1]+z[2]+z[3]+z[4]
    return -4*J*y-h*x

def S(z):
    Sx, Sy, Sz = 0,0,0

    for i in range(6):
        if z[i]<0: return -1
        Sz-=c[i]*xlx(z[i])

    y = [z[0]+2*z[1]+z[3],
        z[1]+z[2]+z[3]+z[4],
        z[2]+2*z[4]+z[5]]

    for i in range(3):
        Sy-=b[i]*xlx(y[i])

    x = [y[0]+y[1], y[1]+y[2]]

    for i in range(2):
        Sx-=a[i]*xlx(x[i])
    
    return Sx-2*Sy+Sz
    
def F(J,h,z,T):
    if T == 0: return H(J,h,z) #might as well skip entropy calc if possible
    return H(J,h,z)-T*S(z)

def min(J, hpj=1):
    h = hpj*abs(J)
    Temp = np.linspace(0,10,201)
    
    #ensure reasonable compositions
    con = opt.LinearConstraint([[ 1, 3, 2, 1, 1, 0],  # 0 <= x   <= 1
                                [ 0, 1, 2, 1, 3, 1],  # 0 <= 1-x <= 1
                                [ 1, 4, 4, 2, 4, 1]], # total proability 1
                               [0,0,1], [1,1,1])
    #restrict 0 <= zi <= 1
    bound = opt.Bounds([0,0,0,0,0,0],[1,1,1,1,1,1])
    
    mF, mX, mY = [],[],[] #I need to choose consistent variable names...
    
    for Tpj in Temp:
        status="Calculating T/J:"+str(Tpj)
        print(status)

        free = lambda z: F(J, h, z, Tpj*abs(J))

        res = opt.minimize(free,                           #minimize free energy
                [0.0625,0.0625, 0.0625, 0.0625, 0.0625, 0.0625], #starting guess at centre
                method = 'trust-constr',                   #optimization method
                constraints = con, bounds = bound)        #set constraints on variables
        
        mF.append(free(res.x))
        mX.append(res.x[0]+3*res.x[1]+2*res.x[2]+res.x[3]+res.x[4])
        mY.append(res.x[0]+2*res.x[1]+res.x[3])

    fig = plt.figure(figsize=plt.figaspect(2.))
    fig.suptitle('2D square lattice, square approx')

    #plot min free energy wrt temp
    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('Temperature (T/J)')
    ax.set_ylabel('Minimum Free Energy')
    ax.plot(Temp, mF)

    #plot composition wrt temp
    ax = fig.add_subplot(2,1,2, projection='3d')
    ax.set_xlabel('Composition: X')
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_zlabel('Temperature (T/J)')
    ax.plot(mX, mY, Temp)

    plt.show()
