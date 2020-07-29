import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import scipy.optimize as opt

#Computations are done in terms of y_ij, which is an asymmetric square matrix; the first index denotes the component on the first (even) sublattice , and the second on the second (odd) sublattice
#the y variable is flattened to form a 1D array to be processed by scipy.optimize

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

def Xo(y):
    xo=[]
    for i in range(3):
        for j in range(3):

    return xo

def Xe(y):
    



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
def min(J=[[]], h=[[]], Trang=[0,4], samp=100):

