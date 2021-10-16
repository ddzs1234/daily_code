#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 13 22:16:47 2018

@author: ashley
"""
from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
#parameter
mu=0.001
m1=0.95459e-3  #m1/ms
m2=0.28581e-3  #m2/ms
a1=5.20336301
a2=9.53707032
e1=0.04839266
e2=0.05415060
w1=14.75385 #degree
w2=92.43194
i1=1.30530 #degree
i2=2.48446
omega1=100.55615 #degree
omega2=113.71504
lamada1=34.40438 #degree
lamada2=49.99432
T1=11.862615
T2=29.447498 

p1=i1*math.sin(omega1*180/math.pi)
p2=i2*math.sin(omega2*180/math.pi)
q1=i1*math.cos(omega1*180/math.pi)
q2=i2*math.cos(omega2*180/math.pi)
h1=e1*math.sin(w1*180/math.pi)
h2=e2*math.cos(w2*180/math.pi)
k1=e1*math.sin(w1*180/math.pi)
k2=e2*math.cos(w2*180/math.pi)


alpha=a1/a2

err=1e-8
tend=6000
#print alpha

def b1(alpha):
    a=alpha
    return 3*a*(1+15*(a**2)/8+175*(a**4)/64+3675*a**(6)/1024+72765*(a**8)/16384)

def b2(alpha):
    a=alpha
    return 15*(a**2)*(1+7*a**2/4+315*a**4/128+1617*a**6/512+63063*(a**8)/16384)/4

def c0(alpha):
    a=alpha
    return 5*a*(1-a**2/8)/4

#print b1(alpha),b2(alpha)

def c10(alpha,b1):
    return 2*math.pi*(m2)*(alpha**2)*b1(alpha)/(4*T1*(1+m1))

def c20(alpha,b1):
    return 2*math.pi*(m1)*(alpha**2)*b1(alpha)/(4*T2*(1+m2))

a11=c10(alpha,b1)*180/math.pi
a12=-c0(alpha)*c10(alpha,b1)*180/math.pi
a21=-c0(alpha)*c20(alpha,b1)*180/math.pi
a22=c20(alpha,b1)*180/math.pi

b11=c10(alpha,b1)*180/math.pi
b12=-c10(alpha,b1)*180/math.pi
b21=-c20(alpha,b1)*180/math.pi
b22=c20(alpha,b1)*180/math.pi

#b1=3.18970
#b2=2.08583
#
#c0=b1/b2
#c1=2*math.pi*(alpha**2)*b1*m2/(4*T1*(1+m1))
#c2=2*math.pi*(alpha**2)*b1*m1/(4*T2*(1+m2))
#
#a11=c1*180/math.pi
#a12=-c0*c1*(180/math.pi)
#a21=-c0*c2*(180/math.pi)
#a22=c2*180/math.pi
#
#b11=c1*180/math.pi
#b12=-c1*180/math.pi
#b21=-c2*180/math.pi
#b22=c2*180/math.pi
#   

b=np.mat([[0,0,0,0,0,0,0,0,0,0,0,0,0],[2/27,0,0,0,0,0,0,0,0,0,0,0,0],\
         [1/36,1/12,0,0,0,0,0,0,0,0,0,0,0],[1/24,0,1/8,0,0,0,0,0,0,0,0,0,0],\
         [5/12,0,-25/16,25/16,0,0,0,0,0,0,0,0,0],[1/20,0,0,1/4,1/5,0,0,0,0,0,0,0,0],\
         [-25/108,0,0,125/108,-65/27,125/54,0,0,0,0,0,0,0],\
         [31/300,0,0,0,61/225,-2/9,13/900,0,0,0,0,0,0],\
         [2,0,0,-53/6,704/45,-107/9,67/90,3,0,0,0,0,0],\
         [-91/108,0,0,23/108,-976/135,311/54,-19/60,17/6,-1/12,0,0,0,0],\
         [2383/4100,0,0,-341/164,4496/1025,-301/82,2133/4100,45/82,45/162,18/41,0,0,0],\
         [3/205,0,0,0,0,-6/41,-3/205,-3/41,3/41,6/41,0,0,0],\
         [-1777/4100,0,0,-341/164,4496/1025,-289/82,2193/4100,51/82,33/164,12/41,0,1,0]])
c=np.mat([41/840,0,0,0,0,34/105,9/35,9/35,9/280,9/280,41/840,0,0])
c1=np.mat([0,0,0,0,0,34/105,9/35,9/35,9/280,9/280,0,41/840,41/840])
a0=np.mat([0,2/27,1/9,1/6,5/12,1/2,5/6,1/6,2/3,1/3,1,0,1])

f=np.mat([[0,0,a11,a12,0,0,0,0],\
         [0,0,a21,a22,0,0,0,0],\
         [a11,a12,0,0,0,0,0,0],\
         [-a21,-a22,0,0,0,0,0,0],\
         [0,0,0,0,0,0,b11,b12],\
         [0,0,0,0,0,0,b21,b22],\
         [0,0,0,0,-b11,-b12,0,0],\
         [0,0,0,0,-b21,-b22,0,0]])

def generate_coef(t,x,h):
    def diff_func(t,x0):
        y=[]
        y=np.dot(f,x0)
        return y
    
    k=np.mat(np.zeros((8,13)))
    for i in range(13):
        k[:,i]=diff_func(t+a0[0,i]*h,x+h*k*b[i].T)
    return k


def RKF78(h,t,x,err,n):
    while True:
        k=generate_coef(t,x,h)
        x0=x+h*k*c.T
        x1=x+h*k*c1.T
        temp=x1-x0
        temp_delta1=(temp[0,0]**2+temp[2,0]**2)**0.5/h
        temp_delta2=(temp[1,0]**2+temp[3,0]**2)**0.5/h
        
        if temp_delta1<err and temp_delta2<err:
            break
        else:
            h=h/2
    return (h,x0)


if __name__=='__main__':
    t=0
    x0=np.mat([-0.011733849701336107,0.039057296306831114,-0.011733849701336107,
               0.039057296306831114,-0.3076450029460795,-0.6658063028887746,
               1.2685277459174105,2.3935838106558083]).T
    Time=[t]
    m1=[[x0[0,0],x0[2,0],x0[4,0],x0[6,0]]]
    m2=[[x0[1,0],x0[3,0],x0[5,0],x0[7,0]]]
    e1=[(x0[0,0]**2+x0[2,0]**2)**0.5]
    e2=[(x0[1,0]**2+x0[3,0]**2)**0.5]
    i1=[(x0[4,0]**2+x0[6,0]**2)**0.5]
    i2=[(x0[5,0]**2+x0[7,0]**2)**0.5]
    while t<tend:
        h=0.01 
        h,x0_new=RKF78(h,t,x0,err,8) 
        x0=x0_new
        t=t+h
        m1.append([x0[0,0],x0[2,0],x0[4,0],x0[6,0]])
        m2.append([x0[1,0],x0[3,0],x0[5,0],x0[7,0]])
        Time.append(t)
        e1.append((x0[0,0]**2+x0[2,0]**2)**0.5)
        e2.append((x0[1,0]**2+x0[3,0]**2)**0.5)
        i1.append((x0[4,0]**2+x0[6,0]**2)**0.5)
        i2.append((x0[5,0]**2+x0[7,0]**2)**0.5)
    m1=np.array(m1)
    m2=np.array(m2)    
    
    plt.figure()
    plt.plot(Time,e1)
    
    plt.figure()
    plt.plot(Time,e2)
    
    plt.figure()
    plt.plot(Time,i1)
    
    plt.figure()
    plt.plot(Time,i2)
        






























    