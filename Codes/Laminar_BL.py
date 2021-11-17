# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 23:24:23 2021

@author: Sidney Williams
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import time

def TDM(n1,T1):
    for ii in range (n1-1):
        j=ii+1
        T1[j,1]=T1[j,1]-(T1[j,0]*T1[j-1,2])/T1[j-1,1]
        T1[j,3]=T1[j,3]-(T1[j,0]*T1[j-1,3])/T1[j-1,1]
        T1[j,0]=0.0
    T1[n1-1,3]=T1[n1-1,3]/T1[n1-1,1]
    for ii in range (n1-1):
        j=(n1-2)-ii
        T1[j,3]=(T1[j,3]-T1[j,2]*T1[j+1,3])/T1[j,1]
        T1[j,1]=1.0 
        T1[j,2]=0.0
    T1[n1-1,1]=1.0
    return(T1)
def maxeps(v,vit,n2,k1,i1):
    maxeps=abs(v[i1,0]-vit[i1,0])
    for ii in range(k1-1):
        j=ii+1
        epsi=abs(v[i1,j]-vit[i1,j])
        if epsi>=maxeps:
            maxeps=epsi
    return(maxeps)

n=2500
k=300

L=0.5
p=3.0
U0=5.0
nu= 1.8e-5
de1=p*L*5.0/math.sqrt(U0*L/nu)
eps=1e-03
M=1.008
b=1.0025
bs=0.0
As=0.0

for i in range(0,n-1):
    bs=bs+b**i
for i in range(0,k-1):
    As=As+M**i
dy1=de1/As
dx1=L/bs

U=np.zeros([n,k], dtype=float)
V=np.zeros([n,k], dtype=float)
Uit=np.zeros([n,k], dtype=float)
Vit=np.zeros([n,k], dtype=float)
abcd=np.zeros([k,4], dtype=float)
x=np.zeros([n], dtype=float)
y=np.zeros([k], dtype=float)
cf1=np.zeros(n, dtype=float)
cf2=np.zeros(n, dtype=float)
Rex=np.zeros(n, dtype=float)
lcf1=np.zeros(n-2, dtype=float)
lcf2=np.zeros(n-2, dtype=float)
lRex=np.zeros(n-2, dtype=float)
U[:,k-1]=U0
U[0,:]=U0
U[:,0]=0.0

for ii in range(0,n-1):
    i=ii+1
    x[i]=x[i-1]+dx1*b**(i-1)
    
for jj in range(0,k-1):
    j=jj+1
    y[j]=y[j-1]+M**(j-1)*dy1

for i in range(1,n):
    U[i,:]=U[i-1,:]
    V[i,:]=U[i-1,:]
    it=0
    dx=x[i]-x[i-1]
    for itt in range(0,400):
        Uit[i,:]=U[i,:]
        Vit[i,:]=V[i,:]
        for jj in range(0,k-1):
            j=jj+1
            dym=y[j]-y[j-1]
            V[i,j]=V[i,j-1]-dym/dx*(U[i,j]-U[i-1,j])
        
        for jj in range(0,k-2):
            j=jj+1
            dym=y[j]-y[j-1]
            dyp=y[j+1]-y[j]
            abcd[j,0]=-V[i,j]*dyp/dym/(dyp+dym)-nu*2.0/(dyp+dym)/dym
            abcd[j,1]=U[i,j]/dx+(dyp-dym)*V[i,j]/dyp/dym+2.0*nu/dyp/dym
            abcd[j,2]=V[i,j]*dym/dyp/(dyp+dym)-nu*2.0/(dyp+dym)/dyp
            abcd[j,3]=U[i,j]*U[i-1,j]/dx
        
        abcd[0,0]=0.0
        abcd[0,1]=1.0
        abcd[0,2]=0.0
        abcd[0,3]=U[i,0]
        
        abcd[k-1,0]=0.0
        abcd[k-1,1]=1.0
        abcd[k-1,2]=0.0
        abcd[k-1,3]=U[i,k-1]
        
        abcd=TDM(k,abcd)
        U[i,:]=abcd[:,3]         
        
        it=itt
        epsitU= maxeps(U,Uit,n,k,i)
        epsitV= maxeps(V,Vit,n,k,i)
        
        if (epsitU <= eps) and (epsitV <= eps):
            break
    Rex[i]=U[i,k-1]*(i+1)*dx/(nu)
    cf1[i]=nu/U[i,k-1]**2.0*(U[i,1]-U[i,0])/(y[1]-y[0])
    cf2[i]=2*0.332/math.sqrt(Rex[i])
for i in range(n-2):
    lRex[i]=math.log(Rex[i+2])
    lcf1[i]=math.log(cf1[i+2])
    lcf2[i]=math.log(cf2[i+2])
#plt.plot(lcf1[:],lRex[:])
#plt.plot(lcf2[:],lRex[:])
plt.plot(y[:],U[n-1,:])
plt.show()
