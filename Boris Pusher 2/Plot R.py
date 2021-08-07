# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 17:47:55 2021

@author: Sidney Williams
"""

import numpy as np
from scipy import special
from scipy import constants
import matplotlib.pyplot as plt

def k(z,r,a):
        kek=np.sqrt(4*r*a/(z**2+(r+a)**2))
        
        return kek
    
def Br(z,r,I,a):
        kek=k(z,r,a)
        mu0=constants.mu_0
        pi=constants.pi
        Br=mu0*I*kek*z/(4*pi*(a*r**3)**(1/2))*(-special.ellipk(kek**2)+(1-1/2*kek**2)/(1+kek**2)*special.ellipe(kek**2))
        
        return Br
    
def Bz(z,r,I,a):
        kek=k(z,r,a)
        mu0=constants.mu_0
        pi=constants.pi
        Bz=mu0*I*kek/(4*pi*(a*r)**(1/2))*(special.ellipk(kek**2)+((a+r)*kek**2-2*r)/(2*r*(1-kek**2))*special.ellipe(kek**2))
        
        return Bz

def B(L,I1,I2,a1,a2,z,r):
        Bzz=Bz(z+L,r,I2,a2)+Bz(z-L,r,I1,a1)
        Brr = Br(z+L,r,I2,a2)+Br(z-L,r,I1,a1)
        B=np.sqrt(Bzz**2+Brr**2)       
            
        return B
    
def R(L,I1,I2,a1,a2,r,r0):
        R=B(L,I1,I2,a1,a2,L,r)/B(L,I1,I2,a1,a2,0,r0)
        
        return R

N=100
r0=0.1
r=0.075
I1=150*1e+06
I2=150*1e+06
a1=0.3
a2=0.2

rangL=np.linspace(0.1,1,N)
rangR=[0]*N
i=0

for L in rangL:
   rangR[i]=R(L,I1,I2,a1,a2,r,r0)
   i += 1
   
plt.plot(rangL,rangR)
plt.show