#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 22:32:14 2022

@author: tongzhu
"""


import numpy as np
from scipy.optimize import brentq 
import matplotlib.pyplot as plt  
from matplotlib.pyplot import figure

import matplotlib.ticker as ticker
import matplotlib as mpl



n = []
xb = []
xs = []
with open("phase_diagram-300K.dat","r") as input_file:
    for line in input_file:
        data = line.split()
        n.append(1.0/float(data[0]))
        n.append(1.0/float(data[0]))
        xb.append(float(data[1]))
        xb.append(float(data[2]))
        xs.append(float(data[3]))
        xs.append(float(data[4]))

x_bin = [x for (x,y) in sorted(zip(xb,n))]
T_bin = [y for (x,y) in sorted(zip(xb,n))]
x_spin = [x for (x,y) in sorted(zip(xs,n))]
T_spin = [y for (x,y) in sorted(zip(xs,n))]



figure(figsize=(6, 6), dpi=600)
plt.tick_params(axis='both', which='both', labelsize=12,direction = 'out') 
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['axes.linewidth'] = 1.8
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['font.family'] = "Arial"
ax = plt.subplot(1,1,1)
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))
ax.tick_params(axis='both',which='both',direction='out',width=1.8,labelsize=22, bottom=True,top=True, left=True,right=True,labelbottom=True,labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis='x',which='both',pad=6.0)
plt.yticks(fontname="Arial")
plt.xticks(fontname="Arial")
plt.yscale("log")


plt.xlabel("Bromine Ratio",fontname= "Arial",fontsize = 22)
plt.ylabel("$\eta$ $\propto$ I$_{sun}$",fontname = "Arial", fontsize =22)





plt.plot(x_bin,T_bin,'b-')
plt.fill_between(x_bin,T_bin,2.0,color='b',alpha=0.3)
plt.plot(x_spin,T_spin,'r-')
plt.fill_between(x_spin,T_spin,2.0,color='r',alpha=0.3)

#font = {'family': 'Times New Roman','size': 28}
plt.xlim([0.0,1.0])
#plt.ylim([0.0,10.0])
plt.ylim([0.05,2.0])
#plt.ylim([200.0,350.0])
plt.xticks(np.arange(0.0,1.1,0.2))
#plt.yticks(np.arange(200.0,375.0,25.0))
#plt.tick_params(axis='both', which='major', labelsize=18)
#plt.xlabel('$x_{\mathregular{Br}}$',fontdict=font)
#plt.ylabel('$T$ (K)',fontdict=font)
plt.savefig('phase_diagram-300K-log.pdf',bbox_inches='tight')
plt.show()
plt.close()