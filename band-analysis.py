#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:31:12 2020

@author: tongzhu
"""

import numpy as np
from scipy.optimize import brentq 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import matplotlib.ticker as ticker
import matplotlib as mpl


# Boltzmann constant in eV per K (CODATA 2006)
kB = 8.617343e-5 

# Variables
# x -> composition of the alloy
# T -> temperature

# Reading input data from file energies.dat
# j -> index of configurations
# nj -> number of atoms B in configuration j
# gj -> degeneracy of configuration j
# ej -> energy of configuration j in eV per formula unit
nj = []
gj = []
ej = []
bv = []
anion_no = 3.0
with open("band-final.dat","r") as input:
    for line in input:
        if line.startswith('#'):
            continue
        data = line.split()
        nj.append(int(data[1]))
        gj.append(int(data[2]))
        bv.append(float(data[4]))
        ej.append(float(data[3])/3.0) # Dividing by 3 in order to obtain energy/anion


# J -> total number of configurations
# N -> number of substitutional atoms in the sublattice
J = len(nj)
N = nj[J-1] 


dbv = []
for j in range(J):
    x = float(nj[j])/float(N)
    dbv.append(bv[j]-(1-x)*bv[0]-x*bv[J-1])





dej = []
X = []
for j in range(J):
    x = float(nj[j])/float(N)
    dej.append(ej[j]-(1-x)*ej[0]-x*ej[J-1])
    X.append(x)

# This quantity is used in the definition of the poliynomium
def alpha(j,x,T):
    return (N*x-nj[j])*gj[j]*np.exp(-dej[j]/(kB*T))

# Poliynomium related to a constraint of GQCA probability
def polynom(eta,x,T):
    total = 0
    for j in range(J):
        total += alpha(j,x,T)*eta**nj[j]
    return total

# Brent method is used to find the root of the polynomial equation
def root(x,T):
    eta_0 = 0.0
    eta_1 = 1e0
    while polynom(eta_1,x,T) >= 0:
        eta_1 = 10*eta_1
    return brentq(polynom,eta_0,eta_1,args=(x,T),maxiter=1000)

# GQCA probability of finding the configuration j in the alloy
def xj(j,T,r0):
    total = 0
    for i in range(J):
        total += gj[i]*r0**nj[i]*np.exp(-dej[i]/(kB*T))
    return gj[j]*r0**nj[j]*np.exp(-dej[j]/(kB*T))/total

# Probability in a random alloy (a priori)
def xj0(j,x):
    return gj[j]*x**nj[j]*(1-x)**(N-nj[j])


def band_V0(x):
    if x == 0 or x == 1:
        return 0
    else:
        total = 0
        for j in range(J):
            total += xj0(j,x)*bv[j]
        return total-(1-x)*bv[0]-x*bv[J-1]



def band_V(x,T):
    if x == 0 or x == 1: 
        return 0
    else:
        r0 = root(x,T)
        total = 0
        for j in range(J):
            total += xj(j,T,r0)*bv[j]
        return total - (1-x)*bv[0]-x*bv[J-1]

# Enthalpy of mixing
def enthalpy(x,T):
    if x == 0 or x == 1:
        return 0
    else:
        r0 = root(x,T)
        total = 0
        for j in range(J):
            total += xj(j,T,r0)*ej[j]
        return total-(1-x)*ej[0]-x*ej[J-1]

# Enthalpy of mixing (random alloy)
def enthalpy0(x):
    if x == 0 or x == 1:
        return 0
    else:
        total = 0
        for j in range(J):
            total += xj0(j,x)*ej[j]
        return total-(1-x)*ej[0]-x*ej[J-1]


# Entropy of mixing
def entropy(x,T):
    if x == 0 or x == 1:
        return 0
    else:
        r0 = root(x,T)
        total = 0
        for j in range(J):
#            print xj0(j,x)
            total += xj(j,T,r0)*np.log(xj(j,T,r0)/xj0(j,x))
        return -kB*(x*np.log(x)+(1-x)*np.log(1-x)+total/N)

# Entropy of mixing (ideal solution)
def entropy0(x):
    if x == 0 or x == 1:
        return 0
    else:
        return -kB*(x*np.log(x)+(1-x)*np.log(1-x))

# Free energy of mixing (Helmholtz)
def free(x,T):
    return enthalpy(x,T)-T*entropy(x,T) 

def d2free(x,T):
    h = 0.00001
    x1 = x-h
    x2 = x+h
    y = free(x,T)
    y1 = free(x1,T)
    y2 = free(x2,T)
    return (y2-2.0*y+y1)/h**2







X2 = np.linspace(0.01,0.99,30)
print d2free(0.5,170)
DIM2 = len(X2)
Y2=np.zeros([DIM2])
Y3=np.zeros([DIM2])
Y4=np.zeros([DIM2])
Y5=np.zeros([DIM2])
Y6=np.zeros([DIM2])
Y7=np.zeros([DIM2])
Y8=np.zeros([DIM2])
Y9=np.zeros([DIM2])
Y10=np.zeros([DIM2])
Y11=np.zeros([DIM2])
for i in range(DIM2):
    Y2[i] = band_V0(X2[i])
    Y3[i] = band_V(X2[i],200)
    Y4[i] = band_V(X2[i],300)
    Y5[i] = band_V(X2[i],300)
    
#    Y5[i] = entropy0(X2[i])*1000
#    Y6[i] = entropy(X2[i],200)*1000
#    Y7[i] = entropy(X2[i],300)*1000
#    Y8[i] = free(X2[i],100)
#    Y9[i] = free(X2[i],200)
#    Y10[i] = free(X2[i],300) 
#    Y11[i] = free(X2[i],400)
#    print Y2[i]


font_label = {'family': 'Arial','size': 25}
font_legend = {'family': 'Arial','size': 18}
font = {'family': 'Arial','size': 28}


figure(figsize=(6, 6), dpi=600)
plt.tick_params(axis='both', which='both', labelsize=12,direction = 'out') 
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['axes.linewidth'] = 1.8
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['font.family'] = "Arial"
ax = plt.subplot(1,1,1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.08))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.016))
ax.tick_params(axis='both',which='both',direction='out',width=1.8,labelsize=22, bottom=True,top=True, left=True,right=True,labelbottom=True,labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis='x',which='both',pad=6.0)
plt.yticks(fontname="Arial")
plt.xticks(fontname="Arial")



plt.xlabel("Bromine Ratio",fontname= "Arial",fontsize = 22)
plt.ylabel("-$\Delta$E$_{light \  contri.}$($x$,T) (eV/formula)",fontname = "Arial", fontsize =22)




#figure(figsize=(6, 4), dpi=600)
#plt.tick_params(axis='both', which='major', labelsize=18,direction = 'in') 
#plt.xlim([0.0,1.0])
#plt.ylim([-0.021,0.003])

plt.plot(X,dbv,'ro')
plt.plot(X2,Y2,'-y')
plt.plot(X2,Y3,'-b')
plt.plot(X2,Y4,'-k')
plt.plot(X2,Y5,'-g')

np.savetxt('band-Contribution-300K.out', (X2,Y5))

#font = {'family': 'Times New Roman','size': 28}
plt.xlim([0.0,1.0])
#plt.axhline(y=0,xmin=0.0,xmax=1.0,color='black')
#plt.ylim([174.0,240.0])
plt.xticks(np.arange(0.0,1.1,0.2))
#plt.yticks(np.arange(174.0,240.0,5.0))
#plt.tick_params(axis='both', which='major', labelsize=18)
#plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
#plt.ylabel('$\Delta E_g $ (eV)',fontdict=font)
plt.savefig('Bandgap.pdf',bbox_inches='tight')
plt.show()
plt.close()

'''

plt.plot(X2,Y5,'-r')
plt.plot(X2,Y7,'-b')
plt.plot(X2,Y6,'-g')

plt.xlim([0.0,1.0])
#plt.axhline(y=0,xmin=0.0,xmax=1.0,color='black')
#plt.ylim([174.0,240.0])
plt.xticks(np.arange(0.0,1.1,0.1))
#plt.yticks(np.arange(174.0,240.0,5.0))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
plt.ylabel('$\Delta S $ (meV k$^{-1}$/anion)',fontdict=font)
plt.savefig('T-S.pdf',bbox_inches='tight')
plt.show()
plt.close()


plt.plot(X2,Y8,'-',color="purple")
plt.plot(X2,Y9,'-b')
plt.plot(X2,Y10,'-g')
plt.plot(X2,Y11,'-',color="grey")

plt.xlim([0.0,1.0])
#plt.axhline(y=0,xmin=0.0,xmax=1.0,color='black')
#plt.ylim([174.0,240.0])
plt.xticks(np.arange(0.0,1.1,0.1))
#plt.yticks(np.arange(174.0,240.0,5.0))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
plt.ylabel('$\Delta F $ (eV/anion)',fontdict=font)
plt.savefig('T-F.pdf',bbox_inches='tight')
plt.show()
plt.close()

#plt.plot(X2,Y8,'-r')
#plt.plot(X2,Y9,'-b')
plt.plot(X2,Y10,'-k')
#plt.plot(X2,Y11,'-k')

plt.ylim([-0.025,0.0])
plt.yticks(np.arange(-0.025,0.001,0.005))


plt.xlim([0.0,1.0])
#plt.axhline(y=0,xmin=0.0,xmax=1.0,color='black')
#plt.ylim([174.0,240.0])
plt.xticks(np.arange(0.0,1.1,0.2))
#plt.yticks(np.arange(174.0,240.0,5.0))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
plt.ylabel('$\Delta F $ (eV/anion)',fontdict=font)
plt.savefig('T-F-300.pdf',bbox_inches='tight')
plt.show()
plt.close()

def eg(x): 
    return 0.38202*x+0.20898*x**2  
def eg2(x):
    return  0.59725103*x*(1-x)




Eg=np.zeros([DIM2])
Final_E = np.zeros([DIM2])
for i in range(DIM2):
    Eg[i] = eg2(X2[i])

nmax = -1.0*Eg/(3.0*Y10)
Final_E = Eg + Y10*39.0

plt.plot(X2,Final_E,'-r')
plt.plot(X2,Y10*39.0,'-k')

plt.xlim([0.0,1.0])
#plt.axhline(y=0,xmin=0.0,xmax=1.0,color='black')
#plt.ylim([-0.025,0.0])
#plt.yticks(np.arange(-0.025,0.001,0.005))
plt.xticks(np.arange(0.0,1.1,0.2))
#plt.yticks(np.arange(174.0,240.0,5.0))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
plt.ylabel('$\Delta F$  (eV)',fontdict=font)
plt.savefig('T-F-photon-n=13.pdf',bbox_inches='tight')
plt.show()
plt.close()



plt.plot(X2,nmax,'r')
plt.axhline(y=0.0,ls='--',color='k')


plt.xlim([0.0,1.0])
#plt.axhline(y=0,xmin=0.0,xmax=1.0,color='black')
#plt.ylim([-0.025,0.0])
#plt.yticks(np.arange(-0.025,0.001,0.005))
plt.xticks(np.arange(0.0,1.1,0.2))
#plt.yticks(np.arange(174.0,240.0,5.0))
plt.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
plt.ylabel('$n_{max}$',fontdict=font)
plt.savefig('domain-size.pdf',bbox_inches='tight')
plt.show()
plt.close()

'''


