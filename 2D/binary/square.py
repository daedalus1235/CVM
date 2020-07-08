import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

#helper functions
def xlx(x):
    if x<0: return 1
    if x==0: return 0
    return x*math.log(x)

def Y(z):
    return [z[0]+2*z[1]+z[2],
            z[1]+z[2]+z[3]+z[4],
            z[2]+2*z[4]+z[5]]
    
def X(z):
    y=Y(z)
    return [y[0]+y[1], y[1]+y[2]]

#Thermodynamic functions
def H(J,h,z):
    x = X(z)

    y = Y(z)
    return -4*J*y[1]-h*x[1]

def S(z):
    #coefficients for multiplicity
    a=[1,1]         #x
    b=[1,2,1]       #y ++ +- --
    c=[1,4,4,2,4,1] #z ++++ +++- ++-- +-+- +--- ----
    
    Sx, Sy, Sz = 0,0,0

    for i in range(len(z)):
        Sz-=c[i]*xlx(z[i])

    y = Y(z)

    for i in range(len(y)):
        Sy-=b[i]*xlx(y[i])

    x = X(z)

    for i in range(len(x)):
        Sx-=a[i]*xlx(x[i])
    
    return Sx-2*Sy+Sz
    
def F(J,h,z,T):
    if T == 0: return H(J,h,z) #might as well skip entropy calc if possible
    return H(J,h,z)-T*S(z)

#minimizatin
def min(J=1, hpj=1):
    h = hpj*abs(J)
    Temp = np.linspace(0,5,51)
    
    #starting guess
    guess=[0.0625,0.0625,0.0625,0.0625,0.0625,0.0625]

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

        res = opt.minimize(free,                   #minimize free energy
                guess,                             #guess is updated to previous result
                method = 'trust-constr',           #optimization method
                constraints = con,                 #set constraints...
                bounds = bound)                    #... and bounds on variables
        
        mF.append(free(res.x))
        tempx=X(res.x)
        mX.append(tempx[0])
        tempy=Y(res.x)
        mY.append(tempy[1])
        guess = res.x
        if tempy[1]>0.5:
            print(res)
        else: print("z: "+str(res.x))
        print("y: "+str(tempy))
        print("x: "+str(tempx))

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
    ax.set_xlim(0,1)
    ax.set_ylabel('Bond Frequency: Y')
    ax.set_ylim(0,0.5)
    ax.set_zlabel('Temperature (T/J)')
    ax.plot(mX, mY, Temp)
    
    plt.show()

