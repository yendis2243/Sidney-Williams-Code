#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 16:51:53 2021

@author: Sidney Williams, divyang
"""

import numpy as np
import pandas as pd
from scipy.optimize import newton
import MainSVW

#%%

class MakeDataFiles:
    
    def __init__(self):
        case="circle"
        
        # Magnetic Mirror, Single Particle 
        if(case=="single"):
         self.x,  self.y,  self.z  = [0], [0], [0]
         self.Vx, self.Vy, self.Vz = [1], [0], [0.5]
         self.M,  self.Q           = [1.672e-25], [1.6e-19]
        
        # Multiple (Identical) Particles In a Circle (const. radius)
        if(case=="circle"):
         Sol = MainSVW.Solver()
         L= 0.1
         I1= 25e+06
         I2= 25e+06
         a1= 0.3
         a2= 0.3
             
         data = pd.read_csv('Inputs.txt', header = None, sep = '=')
         data = data[1].tolist()
         self.N  = int(data[0])
         self.x,  self.y,  self.z  = [0]*self.N, [0]*self.N, [0]*self.N
         self.Vx,  self.Vy,  self.Vz  = [0]*self.N, [0]*self.N, [0.5]*self.N
         self.M,  self.Q  = [1.672e-25]*self.N, [1.6e-19]*self.N
         Vxy=1  #Desired Magnitude of Velocity Perpendicular to the Mag Field
         def initialR(r):
             Br = Sol.Br(L,r,I1,a1)+Sol.Br(-L,r,I2,a2)
             Bz = Sol.Bz(L,r,I1,a1)+Sol.Bz(-L,r,I2,a2)
             root=r*np.sqrt(Br**2+Bz**2)-self.M[0]*Vxy/np.abs(self.Q[0])
             return root
         B0=Sol.Bz(L,0,I1,a1)+Sol.Bz(-L,0,I2,a2)   #Guess B
         x0=0   #Coordinates of the Gyrocenter
         y0=0
         A=0
         for i in range(0,self.N):
             if(A==0):
                 r0=self.M[i]*Vxy/(B0*self.Q[i])
                 A=newton(initialR, r0, fprime=None, args=(), tol=1e-10, maxiter=100, fprime2=None)
                 print(A)
             phi=2*np.pi*i/self.N
             self.x[i]=A*np.sin(phi)+x0
             self.y[i]=A*np.cos(phi)+y0
             self.Vx[i]=-Vxy*self.y[i]/np.sqrt(self.x[i]**2+self.y[i]**2)
             self.Vy[i]=Vxy*self.x[i]/np.sqrt(self.x[i]**2+self.y[i]**2)
        
        # Multiple Particles Energy
        if(case=="energy"):
         minev=1e+03
         maxev=701e+03
         data = pd.read_csv('Inputs.txt', header = None, sep = '=')
         data = data[1].tolist()
         self.N  = int(data[0])
         self.M,  self.Q  = [1.672e-25]*self.N, [1.6e-19]*self.N
         self.x,  self.y,  self.z  = [0]*self.N, [0]*self.N, [0]*self.N
         self.Vx,  self.Vy,  self.Vz  = [0]*self.N, [0]*self.N, [0.5]*self.N
         dev=(maxev-minev)/self.N
         ev=minev
         for i in range(0,self.N):
             KE=ev*1.60218e-19
             self.Vx[i] += np.sqrt(KE*2/self.M[i])
             ev += dev
        
        # Multiple Particles Variable Parallel to Perpendicular Energy Ratio
        if(case=="ratio"):
         minc=0.5
         maxc=4e+06
         KE=1e-14 #Kinetic Energy in Joules
         data = pd.read_csv('Inputs.txt', header = None, sep = '=')
         data = data[1].tolist()
         self.N  = int(data[0])
         self.M,  self.Q  = [1.672e-25]*self.N, [1.6e-19]*self.N
         self.x,  self.y,  self.z  = [0]*self.N, [0.001]*self.N, [0]*self.N
         self.Vx,  self.Vy,  self.Vz  = [0]*self.N, [0]*self.N, [0]*self.N
         dc=(maxc-minc)/self.N
         c=minc
         for i in range(0,self.N):
             self.Vx[i] += np.sqrt(2*KE/(self.M[i]*(1+c)))
             self.Vz[i] += np.sqrt(c)*self.Vx[i]
             c += dc

    def WriteDataFiles(self):
        
        M_data  = open(r'Data/M.txt', 'w')
        Q_data  = open(r'Data/Q.txt', 'w')
        x_data  = open(r'Data/x.txt', 'w')
        y_data  = open(r'Data/y.txt', 'w')
        z_data  = open(r'Data/z.txt', 'w')
        vx_data = open(r'Data/vx.txt','w')
        vy_data = open(r'Data/vy.txt','w')
        vz_data = open(r'Data/vz.txt','w')
        
        np.savetxt(M_data,  np.array([self.M])  )
        np.savetxt(Q_data,  np.array([self.Q])  )
        np.savetxt(x_data,  np.array([self.x])  )
        np.savetxt(y_data,  np.array([self.y])  )
        np.savetxt(z_data,  np.array([self.z])  )
        np.savetxt(vx_data, np.array([self.Vx]) )
        np.savetxt(vy_data, np.array([self.Vy]) )
        np.savetxt(vz_data, np.array([self.Vz]) )
        
        M_data.close()
        Q_data.close()
        x_data.close()
        y_data.close()
        z_data.close()
        vx_data.close()
        vy_data.close()
        vz_data.close()
        

#%%

if __name__ == '__main__':
    
    DoMyWork = MakeDataFiles()
    DoMyWork.WriteDataFiles()
