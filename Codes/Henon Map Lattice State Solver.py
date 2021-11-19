# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 20:52:07 2021

@author: Sidney Williams
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 23:47:22 2021

@author: Sidney Williams
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy.random import rand
from json import dumps
from IPython.display import display, Math
from itertools import product

#Should it produce plots?
plot=True
#Should it be rescaled so a is alone on one side of the equals sign?
Gallas=True
#Size of stretching parameter
a=6
#Maximum length of orbit to be found
n=6
#Length of orbit which will have its Fourier transform taken
n2=6
#Fourier Mode to be displayed
mode=5

#Iterative method for finding the given periodic symbol sequence
def cycle_inverse(symbols):
        signs=[]
        ep=1e-7
        cycle=np.zeros(len(symbols))
        maxiter=1000
        n_iter=0
        deviation=np.zeros(len(symbols))
        deviation[:]=1
        for symbol in symbols:
            if symbol==0:
                signs.append(-1)
            if symbol==1:
                signs.append(1)
        while any(np.greater(abs(deviation),ep)) and n_iter<maxiter:
            n_iter+=1
            for i in range(0,len(symbols)):
                cycle[i]=signs[i]*np.sqrt(abs(1-np.roll(cycle,1)[i]-np.roll(cycle,-1)[i])/a)
            for i in range(0,len(symbols)):
                deviation[i]=np.roll(cycle,-1)[i]-(1-a*(cycle[i])**2-np.roll(cycle,1)[i])
        return(cycle)
        
#Produces the set of all binary sequences up to a specified length which can't
#be built out of shorter binary sequences (01 included, no 0101)
def prime_cycle_itinerary(cyclelength):
        itineraries=[]
        current_cyclelength=1
        for i in range(1,cyclelength+1):
            itinerary=np.zeros(i)
            itinerary[:]=int(1)
            cycle_list=list(product(*((x, 0) for x in itinerary)))
            if itineraries:
                for cycle in cycle_list:
                    add='yes'
                    cyclearray=[]
                    for item in cycle:
                        cyclearray.append(item)
                    if itineraries[0]*current_cyclelength==cycle:
                        add='no'
                    if itineraries[1]*current_cyclelength==cycle:
                        add='no'
                    for item in itineraries:
                        for i in range(0,cyclelength):
                            if list(np.roll(item,i))==cyclearray:
                                add='no'
                            if int(current_cyclelength/len(item))*list(np.roll(item,i))==cyclearray:
                                add='no'
                        if int(current_cyclelength/len(item))*item==cycle:
                            add='no'
                    if add=='yes':
                        itineraries.append(cycle)
            else:
                itineraries.append((0,))
                itineraries.append((1.0,))
            current_cyclelength+=1
        return(itineraries)
        
#Constructs the orbit Jacobian matrix for a specified orbit    
def orbit_jacobian(orbit,orbitlength):
        J=np.arange((orbitlength)**2,dtype='float').reshape(orbitlength,orbitlength)
        c=a
        if(Gallas):
            c=1
        for i in range(0,orbitlength):
            Jrow=[1,2*c*orbit[i],1]
            Jrow.extend(np.zeros(orbitlength-3))
            Ji=np.roll(Jrow,i-1)
            for j in range(0,orbitlength):
                J[i][j]=Ji[j]
        return(J)

#Arranges the specified orbit so that it is centered on its symmetry line if it 
#Has one
def middle_lattice_state(orbit):
        icrit=-1
        icenter=-(len(orbit)//-2)-1
        for i in range(0,len(orbit)):
            if np.abs(orbit[i-1]-orbit[(i+1)%len(orbit)])<=1e-6:
                icrit=i
                break
        
        if icrit==-1:
            for i in range(0,len(orbit)):
                if np.abs(orbit[i]-orbit[(i+1)%len(orbit)])<=1e-6:
                    icrit=i
                    break
        iroll=icenter-icrit
        rolled=np.roll(orbit,iroll)
        return(rolled)
    
       
orbits=[]
fourier=[]
    
scaled=[]
scaled_cycle=[]
fourier_scaled=[]
    
X=[]
Y=[]
 
#Calculates the periodic orbits up to length n   
for cycle in prime_cycle_itinerary(n):
        c=cycle_inverse(cycle)
        for i in range(len(cycle)):
            aprime=np.roll(c,i)
            f=np.fft.fft(aprime)
            fourier.append(f)
        orbits.append(c)

#Calculates the fourier transform    
for cycle in orbits:
        for item in cycle:
            scaled_cycle.append(a*item)
        for i in range(len(cycle)):
            aprime=np.roll(scaled_cycle,i)
            f=np.fft.fft(aprime)
            fourier_scaled.append(f)
        scaled.append(scaled_cycle)
        scaled_cycle=[]

#Rescales which orbits are used    
if(Gallas):
        fouriers=fourier_scaled
        gencycle=scaled
else:
        fouriers=fourier
        gencycle=orbits

#Prepares Fourier modes for plotting    
for cycle in fouriers:
        if len(cycle)==n2:
            X.append(cycle[mode].real)
            Y.append(cycle[mode].imag)

#Calculates the determinant and eigenvalues of the orbit Jacobian for orbits of
#A certain length
dets=[]
eigenvals=[]
cycles_of_interest=[]
eigenvec=[]
for item in gencycle:
        n=len(item)
        if n==5:
            M=orbit_jacobian(middle_lattice_state(item),n)
            eigen=np.linalg.eig(M)
            det=np.linalg.det(M)
            eigenvals.append(eigen[0])
            eigenvec.append(eigen[1])
            cycles_of_interest.append(item)


#Stores all the time reversal symmetric periodic orbits
time_reversible=[]

for cycle in gencycle:
    cyclerev=cycle[::-1]
    for i in range(0,len(cycle)):
        if all(cycle)==all(np.roll(cyclerev,i)):
            time_reversible.append(cycle)
            break

#Plots the given cycle, and its eigenstates, saves the resulting plot
def plotting(cycle):
    n=len(cycle)
    cycle=middle_lattice_state(cycle)
    M=orbit_jacobian(cycle,n)
    eigenvals,eigenvecs=np.linalg.eig(M)
    
    ET=eigenvecs.T
    weights=np.matmul(ET,cycle)
    
    cells=n+1
    l=0
    k=0
    if cells%2==0:
        k=cells/2
        l=2
    else:
        k=cells
        l=1
    fig, axes = plt.subplots(k,l,figsize=(14/2*l,14/3*k))
    
    label=""    
    states=[cycle]
    sites=[]
    colors=[]
    for i in range(0,n):
        if cycle[i]>0:
            label+="1"
        else:
            label+="0"
        if np.abs(cycle[i-1]-cycle[(i+1)%n])<=1e-6:
            colors.append('gold')
        else:
            colors.append('maroon')
        states.append(eigenvecs[:,i])
        sites.append(str(i))
        
    
    iglob=0
    fig.suptitle('Eigenstates of cycle $\overline{'+label+'}$',fontsize=20)
    for i in range(0,k):
        for j in range(0,l):
            if(l==1):
                axes[i].axhline(y=0,color='black')
                axes[i].bar((np.array(sites)),states[iglob],color=colors,width=0.4)
                axes[i].axis(ymin=np.min(states[0]), ymax=np.max(states[0]))
                if iglob==0:
                    axes[i].set_title(r'Cycle: $\overline{'+label+'}$')
                else:
                    axes[i].set_title(r'$\lambda='+str(eigenvals[iglob-1])+'$ \n Weight='+str(weights[iglob-1]))
            else:
                axes[i,j].axhline(y=0,color='black')
                axes[i,j].bar((np.array(sites)),states[iglob],color=colors,width=0.4)
                axes[i,j].axis(ymin=np.min(states[0]), ymax=np.max(states[0]))
                if iglob==0:
                    axes[i,j].set_title(r'Cycle: $\overline{'+label+'}$')
                else:
                    axes[i,j].set_title(r'$\lambda='+str(eigenvals[iglob-1])+'$ \n Weight='+str(weights[iglob-1]))
            
            iglob+=1
    name=str(n)+'latteigvec'+label
    fig.tight_layout(pad=3.0)
    
    plt.savefig(name+'.svg',dpi=300)
    fig.clf()
    
    fig2, axes2 = plt.subplots(1,2,figsize=(14,4.5))
    axes2[0].axhline(y=0,color='black')
    axes2[0].bar((np.array(sites)),weights,color=colors,width=0.4)
    axes2[0].set_title(r'Cycle: $\overline{'+label+'}$ In Eigenbasis')
    
    axes2[1].axhline(y=0,color='black')
    axes2[1].bar((np.array(sites)),states[0],color=colors,width=0.4)
    axes2[1].set_title(r'Cycle: $\overline{'+label+'}$')
    fig2.tight_layout(pad=3.0)
    
    name2=str(n)+'eigenbasis'+label
    plt.savefig(name2+'.svg',dpi=300)
    fig2.clf()

    plt.close("all")
    
    return(weights,eigenvals)
        
    
if(plot):
    eigencoordorbits=np.zeros(len(gencycle)-3)
    eigenvalues=np.zeros(len(gencycle)-3)
    for i in range(1,21):
        eigencoordorbits,eigenvalues=plotting(gencycle[-i])
        print(i)