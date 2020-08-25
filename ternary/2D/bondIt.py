import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import copy

#Computations are done in terms of y_ij, which is an asymmetric square matrix; the first index denotes the component on the first (even) sublattice , and the second on the second (odd) sublattice
#helper functions are written to accept square matrices of any size, so can be used for an n-ary composition of arbitrary n>1.

#helper functions and variables
MAX_ITER    = 2048
Y_PRECISION = 1e-15
T_PRECISION = 1e-9
xA_TOL      = 1e-9

z=4

cA=[[2,1,1],[1,0,0],[1,0,0]]
cB=[[0,1,0],[1,2,1],[0,1,0]]
cC=[[0,0,1],[0,0,1],[1,1,2]]
c=[cA,cB,cC]

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
    return math.sqrt(dot(M,M))

def abstot(M):
    total=0
    for mi in M:
        for mij in mi:
            total+=abs(mij)
    return total            

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
    return add(y, transpose(y))

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
    return z*dot(E,y)

def S(y):
    Sxlx, Syly=0,0
    x=Xe(y)+Xo(y)

    for xi in x:
        Sxlx+=xlx(2*xi)

    for yi in y:
        for yij in yi:
            Syly+=xlx(2*yij)

    return (z-1)/2*Sxlx-z/2*Syly

def F(E,y,T):
    if T==0: return H(E,y)
    return H(E,y)-T*S(y)

#minimization
def ytilde(ycur, Eb, T, m):
    yres=copy.deepcopy(ycur)
    for i in range(len(ycur)):
        for j in range(len(ycur[i])):
            yres[i][j]=math.exp(-Eb[i][j]/T)
            yres[i][j]*=((Xe(ycur)[i]*Xo(ycur)[j])**((z-1)/z))
            #apply lagrange multipliers
            for k in range(len(m)):
                yres[i][j]*=math.exp(m[k]*c[k][i][j]/(z*T))

    yres=normalize(yres,T)  
    return yres

def normalize(ytil, T):
    yres=copy.deepcopy(ytil)
    total = 0
    for yi in ytil:
        for yij in yi:
            total+=yij
    for i in range(len(ytil)):
        for j in range(len(ytil[i])):
            yres[i][j]=ytil[i][j]/(2*total)
    return yres        

def min(Eb=[[0,-1,-1],
            [-1,0,0],
            [-1,0,0]],
        Trang=[0,5],
        samp=200,
        guess=[[25,212.5,212.5],
               [12.5,12.5,0],
               [12.5,0,12.5]],
        m=[1.,0.6,0.5]):
    delta = (Trang[1]-Trang[0])/samp
    temp = np.linspace(Trang[0]+delta, Trang[1], samp)
    tp = np.linspace(Trang[0]+2*delta, Trang[1]-delta, samp-2)
    
    y=normalize(guess,0)

    mY,mF,E,C=[],[],[],[]
    xAe, xBe, xCe =[],[],[]
    xAo, xBo, xCo =[],[],[]
    xAt, xBt, xCt =[],[],[]
    for i in range(len(temp)):
        print('Calculating T='+str(temp[i]))
        counter = MAX_ITER
        T=temp[i]
        yold=[[0, 0, 0],[0, 0, 0],[0, 0, 0]]          #initialize for loop
        change=1

        while change>Y_PRECISION:
            yold=copy.deepcopy(y)
            y=ytilde(yold, Eb, T, m)
            counter-=1
            if(counter<1):
                print('  Max Iterations Reached')
                break
            change=norm(diff(y,yold))

        print('    steps taken:'+str(MAX_ITER-counter))
        print('    last change='+str(change))
        
        x=Xe(y)
        xAe.append(2*x[0])
        xBe.append(2*x[1])
        xCe.append(2*x[2])

        x=Xo(y)
        xAo.append(2*x[0])
        xBo.append(2*x[1])
        xCo.append(2*x[2])
        
        x=Xt(y)
        xAt.append(x[0])
        xBt.append(x[1])
        xCt.append(x[2])

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

    fig = plt.figure(constrained_layout=True)
    fig.suptitle('2D square lattice, bond approx, ternary composition')
    
    gs=fig.add_gridspec(3,2)
    
    ax=fig.add_subplot(gs[0,1])
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Composition')
    ax.set_ylim(-0.05,1.05)
    ax.plot(temp,xAt,label='At', alpha=0.6)
    ax.plot(temp,xBt,label='Bt', alpha=0.6)
    ax.plot(temp,xCt,label='Ct', alpha=0.6)
    ax.legend()

    ax=fig.add_subplot(gs[1,1])
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Free Energy')
    ax.plot(temp,mF)
    
    ax=fig.add_subplot(gs[2,0])
    ax.set_xlabel('Temperature')
    ax.set_ylabel('E=F-T*dF/dT')
    ax.set_ylim(-5,5)
    ax.plot(tp, E)

    ax=fig.add_subplot(gs[2,1])
    ax.set_xlabel('Temperature')
    ax.set_ylabel('C=-T*d2F/dT2')
    ax.set_ylim(-5,5)
    ax.plot(tp, C)
    
    ax=fig.add_subplot(gs[:2,0])
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Composition')
    ax.set_ylim(-0.05,1.05)
    ax.plot(temp,xAe,label='Ae', alpha=0.6)
    ax.plot(temp,xBe,label='Be', alpha=0.6)
    ax.plot(temp,xCe,label='Ce', alpha=0.6)

    ax.plot(temp,xAo,label='Ao', alpha=0.6)
    ax.plot(temp,xBo,label='Bo', alpha=0.6)
    ax.plot(temp,xCo,label='Co', alpha=0.6)
    ax.legend(bbox_to_anchor=(1.05,1), loc='upper left', borderaxespad=0.)

    #high temp slope should be ln(2)~0.693
    slope = mF[samp-1]-mF[math.floor(samp*0.75)]
    slope /= (temp[samp-1]-temp[math.floor(samp*0.75)])
    print('High temp slope: '+str(slope))
    #intercept should be ~0
    intercept=mF[samp-1]-slope*temp[samp-1]
    print('Intercept: ' + str(intercept))

    plt.tight_layout(pad=0, w_pad=0, h_pad=0, rect=(0,0,0.95,0.9))
    plt.show()



#phase diagram
def tsearch(Eb, m, Trang, guess):
    #recursive binary search
    Tavg = (Trang[1]+Trang[0])/2
    
    #base case
    if (Trang[1]-Trang[0])<T_PRECISION:
        return Tavg

    #Find whether to update top or bottom of range
    y=normalize(guess,Tavg)
    yold=[[0,0,0],[0,0,0],[0,0,0]]
    
    counter = MAX_ITER
    change = 1
    exited=False
    
    print('    Tnew='+str(Tavg))
    while change>Y_PRECISION:
        yold=copy.deepcopy(y)
        y=ytilde(yold, Eb, Tavg, m)
        counter-=1
        if(counter<1):
            print('        Max Iterations Reached')
            exited=True
            break
        change=norm(diff(y,yold))
    
    if not exited:
        print('        Iterations taken: '+str(MAX_ITER-counter))

    xdiff=abs(Xe(y)[0]-Xo(y)[0])
    
    if xdiff<=xA_TOL:
        return tsearch(Eb, m, [Trang[0], Tavg], guess)
    return tsearch(Eb, m, [Tavg, Trang[1]], guess)


def phase(Eb=[[0,-1,-1],
            [-1,0,0],
            [-1,0,0]],
        mBC=[1.0,1.0],
        Trang=[0,5],
        mrang=[0,2],
        mnum=200,
        guess=[[25,212.5,212.5],
               [12.5,12.5,0],
               [12.5,0,12.5]],):

    deltamA=mrang[1]-mrang[0]
    deltamA/=mnum
        
    mA_samp=np.linspace(mrang[0]+deltamA,mrang[1],mnum)
    
    Tc=[]

    for mA in mA_samp:
        print('Calculating mA='+str(mA))
        
        m=[mA, mBC[0], mBC[1]]

        Tc.append(tsearch(Eb, m, Trang, guess))

    #output to csv file
    with open('Tc_v_mA_mB'+str(math.floor(10*mBC[0]))+'_mC'+str(math.floor(10*mBC[1]))+'.csv', mode='w') as output:
        outputwriter=csv.writer(output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(mA_samp)):
            outputwriter.writerow([mA_samp[i],Tc[i]])

    plt.plot(mA_samp, Tc)
    plt.title('Transition temperatures with mB='+str(mBC[0])+' and mC='+str(mBC[1]))
    plt.xlabel('mA')
    plt.ylabel('Transition Temperature')

    plt.show()
