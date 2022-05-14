#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:43:07 2022

@author: tongzhu
"""

import numpy as np
from scipy.optimize import brentq 
import matplotlib.pyplot as plt  
from matplotlib.pyplot import figure
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.signal import argrelextrema
from scipy.optimize import fsolve 
Xfree,Free = np.loadtxt('Free-Contribution-350K.out')

Xband,Band = np.loadtxt('band-Contribution-350K.out')


# Cleaning the output file.
open('phase_diagram-350K.dat', 'w').close()


#Initial guess to n=n0+dn
n0 = 0.0
dn = 0.1 
n = n0+dn

#Initial guess to critial compositions 
xc = 0.5 

#Initial guess to the left spinodal point. It will be used to compute the new formula of perovskiters 
xls0 = 0.000001

#Initial guess for binodal points 
#xlb = 0.0000
#xrb = 0.9999

#Headings of the data to be shown 
print " n   left_binodal   right_binodal    left_spinodal   right_spinodal "

t_signal = 1.0
# Search will continue unit the formula variation be small enough
while (dn > 0.01 and t_signal != 0 ): 
    
    #Increase the n 
    n = n0+dn 
    
    #define the free energy(free_light) and its 2nd derivatives(d2free_light) curve according to this n 
    Yfreelight = (3*n*Free-Band)/n
    Fit = np.polyfit(Xfree,Yfreelight,15)
    Fit_curve = np.poly1d(Fit)
    d2_fit = np.polyder(Fit_curve,2)
    
    def free_light(x):
        return Fit_curve(x)
    def d2free_light(x):
        return d2_fit(x)
    
#Testing purpose , checkthe fitting curves         
#    print n
    Xtest = np.linspace(0.0,1.0,100)
    Ytest = free_light(Xtest)
    d2Ytest = d2free_light(Xtest)
    plt.plot(Xtest,Ytest) 
    

# Finding all the spinodal points which d2free should be zero at these points  
# if rootlist beniging/ending equals 0.0/1.0, it stands taht there exist only one spinodal points 
# in the begining and ending regime.
    rootlist=[]
    for i in range(100): 
        root = fsolve(d2free_light,Xtest[i],xtol=1e-12)
 
        if abs(d2free_light(root[0])) < 1e-5:
            if rootlist == []:
                if d2free_light(root[0]+1e-3) > 0 :
                    rootlist.append(0.0)
                    rootlist.append(root[0])
                else:
                    rootlist.append(root[0])
            else:
                signal = 1 
                for j in range(len(rootlist)):
                    if abs(root[0]-rootlist[j]) < 1e-5:
                        signal = 0
                        continue
                if signal == 1:
                    rootlist.append(root[0])
# checing the the last spinodal points d2free is larger than 0 or smaller than 0, 
# if it;s smaller than zero, it means it stands the uphill regime and need another end points for this spinodal region
# which is 1.0 in this case 
    rootlist = sorted(rootlist)
    print n, rootlist
    if rootlist !=[]:               
        if d2free_light(rootlist[-1]+1e-3) < 0:
            rootlist.append(1.0)     
    else:
        #do not find any spinodal points, reset the dn and n values 
        n -= dn/2.0
        break
        
            

#Test purpose, check all the spinodal points.                 
#    print rootlist  
    
    
# pair the whole spinodal points, and then calculate the binodal points for different regions. 
    dimspin = len(rootlist)
    dimpair = dimspin/2  
#    print dimpair      

    binodallist = []        
            
    for i in range(dimpair):
        xls = rootlist[2*i]
        xrs = rootlist[2*i+1]
        
        
        #initial guess of the xc, and binodal points 
        xc = (xls+xrs)/2.0
        if (i==0):
            xlb = 0.0
        else:
            xlb = rootlist[2*i-1]
        
        if (i== dimpair-1):
            xrb = 1.0
        else:
            xrb = rootlist[2*i+2]
        
        #two auxiliary variables 
        x1 = xls
        x2 = xrs 
        
        print xlb,xrb,xls,xrs
        ##Two functions are defined in order to find the common tangent 
        def f1(x):
            return 1e12*((xc-x)*free_light(x2)+(x2-xc)*free_light(x))/(x2-x)
        def f2(x):
            return 1e12*((xc-x1)*free_light(x)+(x-xc)*free_light(x1))/(x-x1)

        # Search for binodal points 
        h1 = 1.0
        h2 = 1.0 
        k  = 0 
        
        if (xls == 0.0): 
            while (h2>1e-6):
                x1 = 0.0
                res = minimize_scalar(f2, bounds=(xrs,xrb), method='bounded')
                h2 = abs(x2-res.x)
                x2 = res.x
                k += 1
                if k == 10:
                    break
        elif (xrs == 1.0):
            while (h1>1e-6):
                x2 = 1.0
                res = minimize_scalar(f1, bounds=(xlb,xls), method='bounded')
                h1 = abs(x1-res.x)
                x1 = res.x
                k += 1
                if k == 10:
                    break
        else:
            while (h1>1e-6 or h2>1e-6):  
                res = minimize_scalar(f1, bounds=(xlb,xls), method='bounded')
 #       print xlb,f1(xlb),free_light(xlb),free_light(x2)
 #       print res.x,f1(res.x),free_light(res.x)
                h1 = abs(x1-res.x)
                x1 = res.x
 #       print xlb,xls,h1,x1
                res = minimize_scalar(f2, bounds=(xrs,xrb), method='bounded')
                h2 = abs(x2-res.x)
#        print h2
                x2 = res.x
                k +=1
                if k == 10:
                    break
        #Binodal points
        xlb = x1
        xrb = x2
        binodallist.append(xlb)
        binodallist.append(xrb)
        
        
        
#    print binodallist
    #showing the results
    #
    print n, rootlist, binodallist
        
    with open("phase_diagram-350K.dat","a") as output:
        output.write("%.2f  " %(n  ))
        for i in range(dimpair):
            output.write("    %.5f %.5f %.5f %.5f    " %(binodallist[2*i],binodallist[2*i+1],rootlist[2*i],rootlist[2*i+1]))
#            output.write("%.2f %.5f %.5f %.5f %.5f\n" %(n,xlb,xrb,xls,xrs))
        output.write("\n" %(  ))       
        #Calculating the new formula n 
    #reset the dn values according to the plots 
    for i in range(dimpair):
        if (rootlist[2*i] != 0.0 and rootlist[2*i+1] !=1.0 ):
            if abs(rootlist[2*i+1]-rootlist[2*i]) < 2e-2: 
                t_signal = 0.0 
                
        #    dn = 0.1*(n-n0)/(xls-xls0)
#    print dn
    n0 = n 
    xls0 = xls 
    
#    dn = 0.1*(n-n0)/(xls-xls0)
    print dn
 
        
#Print the critical points 
        
#print "Critical composition = %.5f" %xc
#print "Critical temperature = %.1f" %n0        
            













