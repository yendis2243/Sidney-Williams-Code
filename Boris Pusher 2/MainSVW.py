# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 16:51:53 2021

@author: Sidney Williams
"""

import pandas as pd
import numpy as np
from scipy import special
from scipy import constants

#%%

class Solver:
    
    def __init__(self):
        
        M    = pd.read_csv(r'Data/M.txt', header = None, sep = ' ')
        Q    = pd.read_csv(r'Data/Q.txt', header = None, sep = ' ')
        x    = pd.read_csv(r'Data/x.txt', header = None, sep = ' ')
        y    = pd.read_csv(r'Data/y.txt', header = None, sep = ' ')
        z    = pd.read_csv(r'Data/z.txt', header = None, sep = ' ')
        vx   = pd.read_csv(r'Data/vx.txt', header = None, sep = ' ')
        vy   = pd.read_csv(r'Data/vy.txt', header = None, sep = ' ')
        vz   = pd.read_csv(r'Data/vz.txt', header = None, sep = ' ')
        data = pd.read_csv('Inputs.txt', header = None, sep = '=')
        data = data[1].tolist()
        
        self.M  =  M.values.tolist()[0]
        self.Q  =  Q.values.tolist()[0]
        self.x  =  x.values.tolist()[0]
        self.y  =  y.values.tolist()[0]
        self.z  =  z.values.tolist()[0]
        self.Vx = vx.values.tolist()[0]
        self.Vy = vy.values.tolist()[0]
        self.Vz = vz.values.tolist()[0]
        self.N  = int(data[0])
        self.T  = float(data[1])
        self.Δt = float(data[2])
        

    def E(self, L, tvary=None):
        
        Ex, Ey, Ez = [0]* self.N, [0]* self.N, [0]* self.N
        
        if(tvary):
            for i in range(0, self.N):
                r = ( self.x[i]**2 + self.y[i]**2 )**0.5
                t = np.arctan2(self.y[i], self.x[i])
            
                Er=-2*r*(L**2+self.z[i]**2)
            
                Ex[i] = Er* np.cos(t)
                Ey[i] = Er* np.sin(t)
                Ez[i] = -2*self.z[i]*r**2

        
        return Ex, Ey, Ez
    
    def F(self,t,tfinal=None,tvary=None):
        Fx, Fy, Fz = [0]* self.N, [0]* self.N, [0]* self.N
        
        # Gravitational Field
        for i in range(0,self.N):
            Fx[i] += 0
            Fy[i] += 0
            g=0
            if(tvary):
                if(t>=tvary):
                    g=9.81*t/tfinal
                    if(t>tfinal):
                        g=9.81
            Fz[i] += -self.M[i]*g
        
        return Fx, Fy, Fz
    
    def k(self,z,r,a):
        kek=np.sqrt(4*r*a/(z**2+(r+a)**2))
        
        return kek
    
    def Br(self,z,r,I,a):
        k=self.k(z,r,a)
        mu0=constants.mu_0
        pi=constants.pi
        Br=mu0*I*k*z/(4*pi*(a*r**3)**(1/2))*(-special.ellipk(k**2)+(1-1/2*k**2)/(1+k**2)*special.ellipe(k**2))
        
        return Br
    
    def Bz(self,z,r,I,a):
        k=self.k(z,r,a)
        mu0=constants.mu_0
        pi=constants.pi
        Bz=mu0*I*k/(4*pi*(a*r)**(1/2))*(special.ellipk(k**2)+((a+r)*k**2-2*r)/(2*r*(1-k**2))*special.ellipe(k**2))
        
        return Bz
    
    def B(self,L,I1,I2,a1,a2):
        
        Bx, By, Bz = [0]* self.N, [0]* self.N, [0]* self.N

        for i in range(0, self.N):
                r = ( self.x[i]**2 + self.y[i]**2 )**0.5
                t = np.arctan2(self.y[i], self.x[i])
            
                Br = self.Br(self.z[i]+L,r,I2,a2)+self.Br(self.z[i]-L,r,I1,a1)
            
                Bx[i] = Br* np.cos(t)
                By[i] = Br* np.sin(t)
                Bz[i] = self.Bz(self.z[i]+L,r,I2,a2)+self.Bz(self.z[i]-L,r,I1,a1)
            
        return Bx, By, Bz
    
    def FirstAdia(self,L,I1,I2,a1,a2):
        Bx, By, Bz = self.B(L,I1,I2,a1,a2)
        B=[0]*self.N
        WPerp=[0]*self.N
        mu=[0]*self.N
        for i in range(0,self.N):
            B[i] += np.sqrt(Bx[i]**2+By[i]**2+Bz[i]**2)
            WPerp[i]=self.M[i]*0.5*(self.Vx[i]**2+self.Vy[i]**2+self.Vz[i]**2-
                 (self.Vx[i]*Bx[i]+self.Vy[i]*By[i]+self.Vz[i]*Bz[i])**2/B[i]**2)
            mu[i] += WPerp[i]/B[i]
        return mu
    
    def Boris(self,L,I1,I2,a1,a2):
        
        ListOfx = open(r'Data/x.txt', 'w')
        ListOfy = open(r'Data/y.txt', 'w')
        ListOfz = open(r'Data/z.txt', 'w')
        
        #ListOfB = open(r'Data/B.txt', 'w')
        #ListOfBr= open(r'Data/Br.txt','w')
        #ListOfBz= open(r'Data/Bz.txt','w')
        
        
        # Step - 1 : Half step velocity
        Ex, Ey, Ez = self.E(L)
        Bx, By, Bz = self.B(L,I1,I2,a1,a2)
        Fx, Fy, Fz = self.F(0)
        
        # Store Previous Positions and Mag Fields
        zp = [0]*self.N            
        diff = [0]*self.N
        
        # First and Second Adiabatic Invariant
        mu0=self.FirstAdia(L,I1,I2,a1,a2)
        J0  = [0]*self.N
        Jf  = [0]*self.N
        Vxp = [0]*self.N
        Vyp = [0]*self.N
        Vzp = [0]*self.N
        B=[0]*self.N
        Br=[0]*self.N
        bounces=[-1]*self.N
        
        for i in range(0, self.N):
            
            QM = self.Q[i]/self.M[i]
            
            self.Vx[i] += (Fx[i]/self.M[i]+QM* (Ex[i] + Bz[i]* self.Vy[i] - By[i]* self.Vz[i]))* self.Δt /2
            self.Vy[i] += (Fy[i]/self.M[i]+QM* (Ey[i] + Bx[i]* self.Vz[i] - Bz[i]* self.Vx[i]))* self.Δt /2
            self.Vz[i] += (Fz[i]/self.M[i]+QM* (Ez[i] + By[i]* self.Vx[i] - Bx[i]* self.Vy[i]))* self.Δt /2
            zp[i] += self.z[i]
            #Vxp[i] += self.Vx[i]
            #Vyp[i] += self.Vy[i]
            #Vzp[i] += self.Vz[i]
        
        
        # Step - 2 : Boris Algorithm
        t  = 0
        qd = [0]* self.N
        hx, hy, hz = [0]* self.N, [0]* self.N, [0]* self.N
        Sx, Sy, Sz = [0]* self.N, [0]* self.N, [0]* self.N
        Ux, Uy, Uz = [0]* self.N, [0]* self.N, [0]* self.N
        Ux_d, Uy_d, Uz_d = [0]* self.N, [0]* self.N, [0]* self.N
        
        for i in range(0, self.N):
            qd[i] = (self.Q[i]* self.Δt)/ (2* self.M[i])
        
        escaped=0
        escape=[0]*self.N
        while (t <= self.T):
            
            for i in range(0, self.N):
                if escape[i]==0:
                    hx[i]   = qd[i]* Bx[i]
                    hy[i]   = qd[i]* By[i]
                    hz[i]   = qd[i]* Bz[i]
                
                    Sx[i]   = 2* hx[i]/ (1 + ( hx[i]**2 + hy[i]**2 + hz[i]**2 ) )
                    Sy[i]   = 2* hy[i]/ (1 + ( hx[i]**2 + hy[i]**2 + hz[i]**2 ) )
                    Sz[i]   = 2* hz[i]/ (1 + ( hx[i]**2 + hy[i]**2 + hz[i]**2 ) )
                
                    Ux[i]   = self.Vx[i] + qd[i]* Ex[i]+self.Δt/(2*self.M[i])*Fx[i]
                    Uy[i]   = self.Vy[i] + qd[i]* Ey[i]+self.Δt/(2*self.M[i])*Fy[i]
                    Uz[i]   = self.Vz[i] + qd[i]* Ez[i]+self.Δt/(2*self.M[i])*Fz[i]
                
                    Ux_d[i] = Ux[i] + Sz[i]* ( Uy[i] + Uz[i]*hx[i] - hz[i]*Ux[i] ) - Sy[i]* ( Uz[i] + Ux[i]*hy[i] - hx[i]*Uy[i] )
                    Uy_d[i] = Uy[i] + Sx[i]* ( Uz[i] + Ux[i]*hy[i] - hx[i]*Uy[i] ) - Sz[i]* ( Ux[i] + Uy[i]*hz[i] - hy[i]*Uz[i] )
                    Uz_d[i] = Uz[i] + Sy[i]* ( Ux[i] + Uy[i]*hz[i] - hy[i]*Uz[i] ) - Sx[i]* ( Uy[i] + Uz[i]*hx[i] - hz[i]*Ux[i] )
                
                    self.Vx[i] = Ux_d[i] + qd[i]* Ex[i]+self.Δt/(2*self.M[i])*Fx[i]
                    self.Vy[i] = Uy_d[i] + qd[i]* Ey[i]+self.Δt/(2*self.M[i])*Fy[i]
                    self.Vz[i] = Uz_d[i] + qd[i]* Ez[i]+self.Δt/(2*self.M[i])*Fz[i]
                
                    self.x[i] += self.Vx[i]* self.Δt
                    self.y[i] += self.Vy[i]* self.Δt
                    self.z[i] += self.Vz[i]* self.Δt
                
                    if (np.sign(self.z[i]-zp[i])!=diff[i]) and (np.abs(self.z[i]-zp[i])>1e-11):
                        bounces[i] += 1
                        print(self.z[i]-zp[i])
                        print(t)
                        wait = input("Press Enter to continue.")
                        diff[i]=np.sign(self.z[i]-zp[i])
                        
                    zp[i]  = self.z[i]
                #Vxp[i] = self.Vx[i]
                #Vyp[i] = self.Vy[i]
                #Vzp[i] = self.Vz[i]
                
                    if (self.z[i]>=L) or (self.z[i]<=-L):
                        escaped += 1
                        escape[i]+=1
                    mu=self.FirstAdia(L,I1,I2,a1,a2)
            
                    if (np.abs(mu0[i]-mu[i])>1e-09):
                        escaped += 1
            
                          
            #np.savetxt(ListOfB, np.array([B]) )
            #np.savetxt(ListOfBr, np.array([Br]) )
            #np.savetxt(ListOfBz, np.array([Bz]) )
                    
            if escaped>=self.N/2:
                print("escaped")
                break
            Ex, Ey, Ez = self.E(L)
            Bx, By, Bz = self.B(L,I1,I2,a1,a2)
            Fx, Fy, Fz = self.F(t)
            
            np.savetxt(ListOfx, np.array([self.x]) )
            np.savetxt(ListOfy, np.array([self.y]) )
            np.savetxt(ListOfz, np.array([self.z]) )
            
            t += self.Δt        
            print(t/self.T)
        ListOfBounces=open(r'Data/bounces.txt','w')
        np.savetxt(ListOfBounces, np.array([bounces]))
        
        ListOfBounces.close()
        ListOfx.close()
        ListOfy.close()
        ListOfz.close()
        
        #ListOfB.close()
        #ListOfBr.close()
        #ListOfBz.close()


#%%

if __name__ == '__main__':
    
    Answer = Solver()
    Answer.Boris(0.1)