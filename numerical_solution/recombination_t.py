#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 13:54:45 2020

@author: ashley
"""
# plot ionization fraction and temperature realtion to calculate recombination temperature, when x~10^-3

import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
from astropy import constants as const



Bh=13.6
yita=10**-10
me=const.m_e
c=const.c
k=const.k_B
e=const.e



t=np.linspace(0.1,1,80)
t=t*e.value/k.value
y=0.26*(5929896575/t)**(3/2)*np.exp(-157821/t)/yita

x=sy.symbols('x')
t_1=[]
x_1=[]
for i in range(0,len(y)):       
    a=sy.solve(x**2/(1-x)-y[i],x)
    for m in a:
        if m>0:
            i_q=t[i]*k.value/e.value
            t_1.append(i_q)
            x_1.append(m)
            

plt.plot(t_1,x_1)
plt.xlabel('T [eV]')
plt.ylabel('ionization fraction x')      
plt.savefig('/home/ashley/Downloads/1.pdf')
