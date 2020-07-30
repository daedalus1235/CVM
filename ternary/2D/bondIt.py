import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

#Computations are done in terms of y_ij, which is an asymmetric square matrix; the first index denotes the component on the first (even) sublattice , and the second on the second (odd) sublattice
#helper functions are written to accept square matrices of any size, so can be used for an n-ary composition of arbitrary n>1.

#helper functions and variables
z=4

def xlx(x):
    if x<0: return 1
    if x==0: return 0
    return x*math.log(x)

def dot(M, N):
    total=0
    for i in range(len(M)):
        for j in range(len(M[i])):
            total+= M[i][j]*N[i][j]
    return total

def norm(M):
    return dot(M,M)

def diff(M,N):
    difference=[]
    for i in range(len(M)):
        di=[]
        for j in range(len(M[i])):
            di.append(M[i][j]-N[i][j])
        difference.append(di)
    return difference

def add(M,N):
    total=[]
    for i in range(len(M)):
        si=[]
        for j in range(len(M[i])):
            si.append(M[i][j]+N[i][j])
        total.append(si)
    return total

def transpose(y):
    yT=[]
    for j in range(len(y[0])):
        yTi=[]
        for i in range(len(y)):
            yTi.append(y[i][j])
        yT.append(yTi)
    return yT

def Yt(y):
    #return add(y, transpose(y))
    return y

def Xe(y):
    xe=[]
    for yi in y:
        xei=0
        for yij in yi:
            xei+=yij
        xe.append(xei)
    return xe    

def Xo(y):
    return Xe(transpose(y))
    
def Xt(y):
    xt=[]
    xe=Xe(y)
    xo=Xo(y)
    for i in range(len(xe)):
        xt.append(xe[i]+xo[i])
    return xt

#thermodynamic functions
def H(E,y):
    return -(z/2)*dot(E,Yt(y))

def S(y):
    Sx,Sy=0,0
    x=Xe(y)#+Xo(y)

    for xi in x:
        Sx-=xlx(xi)

    for yi in y:
        for yij in yi:
            Sy-=xlx(yij)

    return 2*Sy-3*Sx

def F(E,y,T):
    if T==0: return H(E,y)
    return H(E,y)-T*S(y)


#minimization
def min(Eb=[[0,1,1],[1,0,0],[1,0,0]], Trang=[0,4], samp=100):
    T=0
    def iterate():
        yt=y
        for i in range(len(y)):
            for j in range(len(y[i])):
                yt[i][j]=math.exp(-Eb[i][j]/T)
                yt[i][j]*=(Xe(y)[i]*Xo(y)[j])**((z-1)/z)
        
        total = 0
        for yi in yt:
            for yij in yi:
                total+=yij
                
        for i in range(len(yt)):
            for j in range(len(yt[i])):
                y[i][j]=y[i][j]/total
        
        print('Updated y:'+str(yt))

        return yt
    
    delta = (Trang[1]-Trang[0])/samp
    temp = np.linspace(Trang[0]+delta, Trang[1], samp)
    tp = np.linspace(Trang[0]+2*delta, Trang[1]-delta, samp-2)

    mY,mF,E,C=[],[],[],[]
    
    for i in range(len(temp)):
        print('Calculating T='+str(temp[i]))
        counter = 100
        T=temp[i]
        y=[[0.333, 0, 0],[0, 0.333, 0],[0, 0, 0.333]] #starting guess
        yt=[[0, 0, 0],[0, 0, 0],[0, 0, 0]]            #initialize

        while(norm(diff(y,yt))>1e-9):
            y=yt
            counter-=1
            if(counter<1):
                break
            yt=iterate()
        
        mY.append(y)
        mF.append(F(Eb,y,T))

        if i>1:
            #calculate E
            dF=mF[i]-mF[i-2]
            dF/=(delta*2)
            E.append(mF[i-1]-temp[i-1]*dF)

            #calculate C
            ddF=mF[i]-2*mF[i-1]+mF[i-2]
            ddF/=(delta**2)
            C.append(-temp[i-1]*ddF)

    fig = plt.figure()
    fig.suptitle('2D square lattice, bond approx, ternary composition')
    
    #todo: composition plot

    ax=fig.add_subplot(2,2,2)
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Free Energy')
    ax.plot(temp,mF)

    ax=fig.add_subplot(2,2,3)
    ax.set_xlabel('Temperature')
    ax.set_ylabel('E=F-T*dF/dT')
    ax.set_ylim(-5,5)
    ax.plot(tp, E)

    ax=fig.add_subplot(2,2,4)
    ax.set_xlabel('Temperature')
    ax.set_ylabel('C=-T*d2F/dT2')
    ax.set_ylim(-5,5)
    ax.plot(tp, C)

    #high temp slope should be ln(2)~0.693
    slope = mF[samp-1]-mF[math.floor(samp*0.75)]
    slope /= (temp[samp-1]-temp[math.floor(samp*0.75)])
    print('High temp slope: '+str(slope))
    #intercept should be ~0
    intercept=mF[samp-1]-slope*temp[samp-1]
    print('Intercept: ' + str(intercept))

    plt.show()






