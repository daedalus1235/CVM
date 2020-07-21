import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

#variables
A=0
B=1
C=2

#helper functions
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

def Xa(y, comp):
    total=0
    for i in range(int(math.sqrt(len(y)))):
        total+=y[3*comp + i]
    return total

def Xb(y, comp):
    total=0
    for i in range(int(math.sqrt(len(y)))):
        total+=y[comp+3*i]
    return total

def X(y):
    x=[]
    xi=[]
    for i in range(int(math.sqrt(len(y)))):
        xi.append(Xa(y, i))
    x.append(xi)
    xi=[]
    for i in range(int(math.sqrt(len(y)))):
        xi.append(Xb(y, i))
    x.append(xi)
    return x

def tot(y):
    x=X(y)
    total = 0
    for xi in x: 
        for xij in xi:
            total += xij
    return total

#thermodynamic functions
def H(J,h,y):
    return -dot(J,y)-dot(h,X(y))

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

def F(J,h,y,T):
    if T==0: return H(J,h,y)
    return H(J,h,y)-T*S(y)

#minimization
def min(J=[[1,1,1],
            [1,1,1],
            [1,1,1]],
        h=[[0,1,1 ],
            [1,0,0 ]],
        samp=51,
        Trang=[0,4]):

    Temp = np.linspace(Trang[0], Trang[1], samp+1)

    #starting guess
    guess=[]
    for i in range(9):
        guess.append(1/18)

    con = opt.LinearConstraint([[ 2, 2, 2, 2, 2, 2, 2, 2, 2],  #total probability 1
                                [ 2, 2, 2, 0, 0, 0, 0, 0, 0],  #0<x1a<0.5
                                [ 0, 0, 0, 2, 2, 2, 0, 0, 0],  #0<x1b<0.5
                                [ 0, 0, 0, 0, 0, 0, 2, 2, 2],  #0<x1c<0.5
                                [ 2, 0, 0, 2, 0, 0, 2, 0, 0],  #0<x2a<0.5
                                [ 0, 2, 0, 0, 2, 0, 0, 2, 0],  #0<x2b<0.5
                                [ 0, 0, 2, 0, 0, 2, 0, 0, 2]]  #0<x2c<0.5
                               [1,0,0,0,0,0,0],[1,0.5,0.5,0.5,0.5,0.5,0.5])

    bound = opt.bounds([0,0,0,0,0,0,0,0,0],[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])

    mF, E, C =[],[],[]
    
