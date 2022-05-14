#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:18:12 2020

@author: tongzhu
"""

import numpy as np
from scipy.optimize import brentq 
import matplotlib.pyplot as plt  
from matplotlib.pyplot import figure

import matplotlib.ticker as ticker
import matplotlib as mpl

Xfree,Free = np.loadtxt('Free-Contribution-300K.out')

Xband,Band = np.loadtxt('band-Contribution-300K.out')


n = 9



Free_total = 3*n*Free 

Total = Free_total-Band

#font = {'family': 'Times New Roman','size': 28}
font_label = {'family': 'Arial','size': 25}
font_legend = {'family': 'Arial','size': 18}
font = {'family': 'Arial','size': 28}
#figure(figsize=(6, 4), dpi=600)
#plt.tick_params(axis='both', which='major', labelsize=18,direction = 'in') 

figure(figsize=(6, 6), dpi=600)
plt.tick_params(axis='both', which='both', labelsize=12,direction = 'out') 
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['axes.linewidth'] = 1.8
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['font.family'] = "Arial"
ax = plt.subplot(1,1,1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.01))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.002))
ax.tick_params(axis='both',which='both',direction='out',width=1.8,labelsize=22, bottom=True,top=True, left=True,right=True,labelbottom=True,labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis='x',which='both',pad=6.0)
plt.yticks(fontname="Arial")
plt.xticks(fontname="Arial")



plt.xlabel("Bromine Ratio",fontname= "Arial",fontsize = 22)
plt.ylabel("$\Delta$F($\eta,x$,T) (eV/formula)",fontname = "Arial", fontsize =22)

#plt.text()










plt.xlim([0.0,1.0])
plt.xticks(np.arange(0.0,1.1,0.2))

ystring = '$\Delta F(n,x) (eV)$ '
Name = 'T-F-photon-n='+str(n)

labelstring = 'n='+str(n)

#plt.ylim([-0.01,0.02])
#plt.yticks(np.arange(-0.01,0.021,0.01))
#plt.xticks(np.arange(0.0,1.1,0.2))
#plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
#plt.ylabel(ystring,fontdict=font)
plt.tick_params(axis='both', which='major', labelsize=18)

#plt.title(labelstring,fontdict=font)

plt.plot(Xfree,Free_total/n,'k-',label='Dark') 
plt.plot(Xfree,Total/n,'r-')
plt.savefig(Name,bbox_inches='tight')


'''
#begin to get the binodal and spinodal points 
X=Xfree 
Y=Total/n 


Fit=np.polyfit(X,Y,15)
print Fit
Fit_curve= np.poly1d(Fit)
X_fit = np.linspace(0.0,1.0,400)
plt.plot(X_fit,Fit_curve(X_fit),'b--')

Deriv = np.polyder(Fit_curve)
two_Deriv = np.polyder(Fit_curve,2)


X_test = np.linspace(0.0,1.0,10000) 

for i in X_test:
    if abs(two_Deriv(i)) < 1.0e-6:
        print two_Deriv(i)
    else:
        print "No spinodal points"
        


#print(Deriv(X_fit))


'''




#plt.plot(xfree, )