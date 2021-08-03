# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 16:51:53 2021

@author: Sidney Williams
"""

import pandas as pd
import numpy as np

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
        
        # Case - 1 : Constant field
        # for i in range(0, self.N):
        #     Ex[i] += 1
        #     Ey[i] += 0
        #     Ez[i] += 0
        
        # Case - 2 : Nonuniform E Field
        #for i in range(0, self.N):
        #    Ex[i] += 1* np.cos(2* self.x[i])
        #    Ey[i] += 0
        #    Ez[i] += 0
            
        # Case - 3 : Interaction field
        # for i in range(0, self.N):
        #     for j in range(0, self.N):
        #         if (i != j):
                    
        #             r3 = ( (self.x[i] - self.x[j])**2 + (self.y[i] - self.y[j])**2 + (self.z[i] - self.z[j])**2 )**1.5
                    
        #             Ex[i] += self.Q[j]* (self.x[i] - self.x[j])/ r3
        #             Ey[i] += self.Q[j]* (self.y[i] - self.y[j])/ r3
        #             Ez[i] += self.Q[j]* (self.z[i] - self.z[j])/ r3
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
    
    def B(self,L):
        
        Bx, By, Bz = [0]* self.N, [0]* self.N, [0]* self.N
        
        # Case - 1 : Constant B field
        # Bz = [1]* self.N
        
        # Case - 2 : Constant Grad - B
        # for i in range(0, self.N):
        #     Bz[i] = 3 + self.y[i]* 2
        
        # Case - 3 : Constant B theta field
        # Bt = 3
        
        # for i in range(0, self.N):
        #     t = np.arctan2(self.y[i], self.x[i])
            
        #     Bx[i] = - Bt* np.sin(t)
        #     By[i] =   Bt* np.cos(t)
        
        # Case - 4 : Magnetic Mirror
        
        Bzo    = 1

        for i in range(0, self.N):
            #if self.z[i]>=0:
                r = ( self.x[i]**2 + self.y[i]**2 )**0.5
                t = np.arctan2(self.y[i], self.x[i])
            
                Br = -Bzo*r*self.z[i]/L**2
            
                Bx[i] = Br* np.cos(t)
                By[i] = Br* np.sin(t)
                Bz[i] = Bzo*(1 + self.z[i]**2/L**2)
            
            #if self.z[i]<0:
            #    r = ( self.x[i]**2 + self.y[i]**2 )**0.5
            #    t = np.arctan2(self.y[i], self.x[i])
            #
            #    Br = Bzo*r/(2*L)
            #
            #    Bx[i] = Br* np.cos(t)
            #    By[i] = Br* np.sin(t)
            #    Bz[i] = Bzo*(1 - self.z[i]/L)
             
        return Bx, By, Bz
    
    def FirstAdia(self,L):
        Bx, By, Bz = self.B(L)
        B=[0]*self.N
        WPerp=[0]*self.N
        mu=[0]*self.N
        for i in range(0,self.N):
            B[i] += np.sqrt(Bx[i]**2+By[i]**2+Bz[i]**2)
            WPerp[i]=self.M[i]*0.5*(self.Vx[i]**2+self.Vy[i]**2+self.Vz[i]**2-
                 (self.Vx[i]*Bx[i]+self.Vy[i]*By[i]+self.Vz[i]*Bz[i])**2/B[i]**2)
            mu[i] += WPerp[i]/B[i]
        return mu
    
    def SecondAdia(self,L,i,Bxp,Byp,Bzp,Vxp,Vyp,Vzp,zp):
        Bx, By, Bz=self.B(L)
        dx=zp-self.z[i]
        J = dx*self.M[i]/2*((Vxp*Bxp+Vyp*Byp+Vzp*Bzp)/
            np.sqrt(Bxp**2+Byp**2+Bzp**2)+(self.Vx[i]*Bx[i]+self.Vy[i]*By[i]+self.Vz[i]*Bz[i])/
            np.sqrt(Bx[i]**2+By[i]**2+Bz[i]**2))
        
        return J
    
    def Boris(self,L):
        
        ListOfx = open(r'Data/x.txt', 'w')
        ListOfy = open(r'Data/y.txt', 'w')
        ListOfz = open(r'Data/z.txt', 'w')
        
        #ListOfB = open(r'Data/B.txt', 'w')
        #ListOfBr= open(r'Data/Br.txt','w')
        #ListOfBz= open(r'Data/Bz.txt','w')
        
        
        # Step - 1 : Half step velocity
        Ex, Ey, Ez = self.E(L)
        Bx, By, Bz = self.B(L)
        Fx, Fy, Fz = self.F(0)
        
        # Store Previous Positions and Mag Fields
        Bxp, Byp, Bzp = self.B(L)
        zp = [0]*self.N            
        diff = [0]*self.N
        
        # First and Second Adiabatic Invariant
        mu0=self.FirstAdia(L)
        J0  = [0]*self.N
        Jf  = [0]*self.N
        Vxp = [0]*self.N
        Vyp = [0]*self.N
        Vzp = [0]*self.N
        B=[0]*self.N
        Br=[0]*self.N
        
        for i in range(0, self.N):
            
            QM = self.Q[i]/self.M[i]
            
            self.Vx[i] += (Fx[i]/self.M[i]+QM* (Ex[i] + Bz[i]* self.Vy[i] - By[i]* self.Vz[i]))* self.Δt /2
            self.Vy[i] += (Fy[i]/self.M[i]+QM* (Ey[i] + Bx[i]* self.Vz[i] - Bz[i]* self.Vx[i]))* self.Δt /2
            self.Vz[i] += (Fz[i]/self.M[i]+QM* (Ez[i] + By[i]* self.Vx[i] - Bx[i]* self.Vy[i]))* self.Δt /2
            zp[i] += self.z[i]
            Vxp[i] += self.Vx[i]
            Vyp[i] += self.Vy[i]
            Vzp[i] += self.Vz[i]
        
        
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
        count=0
        changes = [0]*self.N
        while (t <= self.T):
            
            for i in range(0, self.N):
            
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
                
                #if (np.sign(self.z[i]-zp[i])!=diff[i]):
                #    changes[i] += 1
                #    print(changes[i])
                #    wait = input("Press Enter to continue.")
                #    diff[i]=np.sign(self.z[i]-zp[i])
                
                #if changes==2:
                #    J0[i] += self.SecondAdia(self,L,i,Bxp[i],Byp[i],Bzp[i],Vxp[i],Vyp[i],Vzp[i],zp[i])
                #elif changes==3 and count==0:
                #    count += 1
                #    print(J0[i])
                #    wait = input("Press Enter to continue.")
                    
                #zp[i]  = self.z[i]
                #Vxp[i] = self.Vx[i]
                #Vyp[i] = self.Vy[i]
                #Vzp[i] = self.Vz[i]
                
                #B[i]=np.sqrt(Bx[i]**2+By[i]**2+Bz[i]**2)
                #Br[i]=np.sqrt(Bx[i]**2+By[i]**2)
                
                if (self.z[i]>=L) or (self.z[i]<=-L):
                    escaped += 1
                mu=self.FirstAdia(L)
            
                if (np.abs(mu0[i]-mu[i])>1e-09):
                    escaped += 1
            
            #np.savetxt(ListOfB, np.array([B]) )
            #np.savetxt(ListOfBr, np.array([Br]) )
            #np.savetxt(ListOfBz, np.array([Bz]) )
                    
            if escaped>=self.N/2:
                print("escaped")
                break
            Ex, Ey, Ez = self.E(L)
            Bx, By, Bz = self.B(L)
            Fx, Fy, Fz = self.F(t)
            
            np.savetxt(ListOfx, np.array([self.x]) )
            np.savetxt(ListOfy, np.array([self.y]) )
            np.savetxt(ListOfz, np.array([self.z]) )
            
            t += self.Δt        
            print(t/self.T)
            
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