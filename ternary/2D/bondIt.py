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

def norm(y):
    return dot(y,y)

def transpose(y):
    yT=[]
    for j in range(len(y[0])):
        yTi=[]
        for i in range(len(y)):
            yTi.append(y[i][j])
        yT.append(yTi)
    return yT

def Yt(y):
    yT=transpose(y)
    Yt=[]
    for i in range(len(y)):
        Yti=[]
        for j in range(len(y[i])):
            Yti.append(y[i][j]+y[j][i])
        Yt.append(Yti)    
    return Yt        

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
def H(J,h,y):
    return -(z/2)dot(E,Yt(y))

def S(y):
    Sx,Sy=0,0
    x=X(y)

    for xi in x:
        for xij in xi:
            Sx-=xlx(xij)

    for yi in y:
        for yij in yi:
            Sy-=xlx(yij)

    return 2*Sy-3*Sx

def F(E,y,T):
    if T==0: return H(E,y)
    return H(E,y)-T*S(y)


#minimization
def min(E=[[]], Trang=[0,4], samp=100):
    T=0

    def iterate():
        yt=y
        for i in range(len(y)):
            for j in range(len(y[i])):
                yt[i][j]=math.exp(-E[i][j]/T) * (Xe(y)*Xo(y))**((z-1)/z)

    
    
    temp = np.linspace(Trang[0], Trang[1], samp+1)
    for i in range(len(temp)):
        counter = 100
        T=
        while(norm(yT)>1):
            iterate()
            counter-=1
            if(counter<1):
                break

    return 0
