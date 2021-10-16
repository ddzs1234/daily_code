#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 16:53:32 2018

@author: ashley
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#import mpl_toolkits.mplot3d.axes3d as p3

mu=0.1
x=np.linspace(-1.5,1.5,45)
X,Y=np.meshgrid(x,x)


def f(x,y):
    r1=((x+mu)**2+y**2)**0.5
    r2=((x+mu-1)**2+y**2)**0.5
    omega=0.5*(x**2+y**2)+(1-mu)/r1+mu/r2
    return 2*omega



c=plt.contour(X, Y, f(X, Y),400,cmap='RdGy')
#print(f(X,Y))
#plt.imshow(f(X,Y),cmap='RdGy')
plt.colorbar()
plt.axis(aspect='image')
#plt.contourf(X, Y, f(X, Y),alpha=0.5,colors='r',cmap=plt.cm.hot)


#plt.xlim(-2,2)  
#plt.ylim(-2,2) 






















