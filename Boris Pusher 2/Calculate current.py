# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 16:13:32 2021

@author: Sidney Williams
"""

import pandas as pd
import numpy as np
from scipy import special
from scipy import constants

def k(z,r,a):
        kek=np.sqrt(4*r*a/(z**2+(r+a)**2))
        
        return kek
    
def Bz(z,r,a):
        kek=k(z,r,a)
        mu0=constants.mu_0
        pi=constants.pi
        Bz=mu0*kek/(4*pi*(a*r)**(1/2))*(special.ellipk(kek**2)+((a+r)*kek**2-2*r)/(2*r*(1-kek**2))*special.ellipe(kek**2))
        
        return Bz
 
m=1.672e-25
q=1.6e-19
J=9.6131e-14
t=0.95
L=1

I=m*np.sqrt(2*J/(m*(1+t)))/(0.075*q*(Bz(L,0.1,0.3)+Bz(L,0.1,0.2)))   