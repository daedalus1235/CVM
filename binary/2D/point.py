import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

Z=4

def xlx(x):
    if x<0:return 1
    if x==0:return 0
    if x>0.5:
        return x*math.log1p(x-1)
    return x*math.log(x)

def H(J,h,x):
    return -J*(Z/2)*(4*x*(1-x)-1)-h*x

def S(x):
    return -xlx(x)-xlx(1-x)

def F(J,h,x,T):
    if T==0: return H(J,h,x)
    return H(J, h, x)-T*S(x)

def min(J=-1, hpj=0, samp=200, Trang=[0,5]):
    Temp = np.linspace(Trang[0],Trang[1],samp+1) 
    
    bound = opt.Bounds([0.0],[1.0])
    guess = [0.]
    
    Free, E, C, mX=[],[],[],[]
    delta = (Trang[1]-Trang[0])/samp
    for i in range(len(Temp)):
        working="Calculating: T/J="+str(Temp[i])
        print(working)
        fun = lambda x: F(J, hpj*abs(J), x, Temp[i]*abs(J))
        res=opt.minimize(fun, guess, method='trust-constr', bounds=bound)
        mX.append(res.x)
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
    fig.suptitle('2D square lattice, point approx')
        
    #plot min free energy wrt temp
    print('Plotting free energy')
    ax = fig.add_subplot(2,2,1)
    ax.set_xlabel('Temperature (T/J)')
    ax.set_ylabel('Min Free Energy')
    ax.plot(Temp, Free, label='min free energy')
   
 
    #plot composition wrt temp
    print('Plotting composition')
    ax=fig.add_subplot(2,2,2)
    ax.plot(mX, Temp, label='min')
    ax.set_xlabel('Composition: X')
    ax.set_xlim(0.0,1.0)
    ax.set_ylabel('Temperature (T/J)')
    ax.set_ylim(0,Trang[1])

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

    with open('point.csv', mode='w') as output:
        outputwriter = csv.writer(output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(Temp)):
            outputwriter.writerow([Temp[i], mX[i]])

    #display plots
    plt.show()
