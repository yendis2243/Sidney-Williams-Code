#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.random import normal
import pandas as pd

#%%

class MakeDataFiles:
    
    def __init__(self):
        
         # Example - 7 : Magnetic Mirror (Nonuniform B field)         
         #self.x,  self.y,  self.z  = [0], [0], [0]
         #self.Vx, self.Vy, self.Vz = [0], [1], [0.5]
         #self.M,  self.Q           = [1.672e-25], [1.6e-19]
        
        # Multiple Particles In a Circle
        # data = pd.read_csv('Inputs.txt', header = None, sep = '=')
        # data = data[1].tolist()
        # self.N  = int(data[0])
        # self.x,  self.y,  self.z  = [0]*self.N, [0]*self.N, [0]*self.N
        # self.Vx,  self.Vy,  self.Vz  = [0]*self.N, [0]*self.N, [0.5]*self.N
        # self.M,  self.Q  = [1.672e-25]*self.N, [1.6e-19]*self.N
        # Vxy=1
        # B0=1
        # x0=0
        # y0=0
        # for i in range(0,self.N):
        #     A=self.M[i]*Vxy/(B0*self.Q[i])
        #     phi=2*np.pi*i/self.N
        #     self.x[i]=A*np.sin(phi)+x0
        #     self.y[i]=A*np.cos(phi)+y0
        #     self.Vx[i]=-Vxy*self.y[i]/np.sqrt(self.x[i]**2+self.y[i]**2)
        #     self.Vy[i]=Vxy*self.x[i]/np.sqrt(self.x[i]**2+self.y[i]**2)
        
        # Multiple Particles Energy
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
