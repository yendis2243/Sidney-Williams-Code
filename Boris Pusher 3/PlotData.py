#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 17:01:44 2020

@author: divyang
"""

from numpy import array, transpose, sum, append, delete, shape, sqrt, linspace, sin
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D

#%%

class Energy:

    '''
        This class will calculate Kinetic, Potential, and Total Energy per particle.
    '''

    def __init__(self, Q, M, N, x, y, z, Vx, Vy, Vz):

        self.Q = Q
        self.M,  self.N  = M,  N
        self.x,  self.y,  self.z  = x,  y,  z
        self.Vx, self.Vy, self.Vz = Vx, Vy, Vz


    def Pot(self, x, y, z):
        
        W = 0

        #for i in range(0, self.N):
        
            # Case - 1 : External Force Potential
            # W -= x[i]* 1 + y[i]* 0 + z[i]* 0
            
            # Case - 2 : Nonuniform E field
            #W -= 0.5* sin(2* x[i])
            
            # Case - 3 : Interaction Potential
            # for j in range(0, self.N):
            #     if (i != j):

            #         r = ( (x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2 )**0.5

            #         W += 0.5* self.Q[i]* self.Q[j]/ r

        return W


    def EnergyCalculations(self):

        # kinetic energy
        KE = 0
        for i in range(0, len(self.M)):
            KE += 0.5* self.M[i]* ( self.Vx[:,i]**2 + self.Vy[:,i]**2 + self.Vz[:,i]**2 )

        # electrostatic potential energy
        p, _ = shape( self.x )
        PE = array([])
        for i in range(0, p):
            PE_dt = self.Pot( self.x[i,:], self.y[i,:], self.z[i,:] )
            PE = append(PE, PE_dt)

        # Total energy
        E = KE + PE

        # Energy per particle
        E  = E / self.N
        KE = KE/ self.N
        PE = PE/ self.N

        return E, KE, PE

#%%

class PlotData:

    def __init__(self):

        # Read Data Files
        M    = pd.read_csv(r'Data/M.txt',  header = None, sep = ' ')
        Q    = pd.read_csv(r'Data/Q.txt',  header = None, sep = ' ')
        x    = pd.read_csv(r'Data/x.txt',  header = None, sep = ' ')
        y    = pd.read_csv(r'Data/y.txt',  header = None, sep = ' ')
        z    = pd.read_csv(r'Data/z.txt',  header = None, sep = ' ')
        data = pd.read_csv('Inputs.txt', header = None, sep = '=')
        data = data[1].tolist()

        self.M  = M.values.tolist()[0]
        self.Q  = Q.values.tolist()[0]
        self.N  = int(data[0])
        self.Δt = float(data[2])
        self.x  = array(x)
        self.y  = array(y)
        self.z  = array(z)
        
        p, _ = shape(self.x)
        self.time = linspace(0, (p-1)* data[2], p)

    def AnimamteParticles(self, Skip = 1):

        gs = GridSpec(1, 2)

        fig = plt.figure( figsize = (80, 60) )
        ax1 = fig.add_subplot( gs[0, 0], projection = '3d' )
        ax2 = fig.add_subplot( gs[0, 1] )
        
        zs=[]
        ts=[]

        for frames, t in enumerate( self.time[::Skip] ):

            fS = frames* Skip

            # PLOT - 1: 3D Particle Position
            #ax1.clear()

            ax1.scatter( self.x[fS], self.y[fS], self.z[fS], marker = 'o' )

            ax1.set_xlabel('x', fontsize = 16)
            ax1.set_ylabel('y', fontsize = 16)
            ax1.set_zlabel('z', fontsize = 16)
            ax1.set_xlim([-0.1, 0.1])
            ax1.set_ylim([-0.1, 0.1])
            ax1.set_zlim([-0.1, 0.1])

            # PLOT - 2: 2D Particle Postion
            
            zs.append(self.z[fS])
            ts.append(t)

            ax2.plot(ts, zs,color="black")
            ax2.set_xlabel('t', fontsize = 16)
            ax2.set_ylabel('z', fontsize = 16)

            plt.pause(1e-11)



#%%

if __name__ == '__main__':

    DoMyWork = PlotData()
    DoMyWork.AnimamteParticles()
