# -*- coding: utf-8 -*-
"""
Created on Thu May 17 18:35:17 2018

@author: surface
"""

'''用自适应积分方法计算该积分，从而得到laplace系数'''
import math
import numpy as np
import matplotlib.pyplot as plt
err=1e-10    #此处为计算积分总的误差

def laplace(j,s,alpha,err):
    f1=lambda x:math.cos(j*x)/(1-2*alpha*math.cos(x)+alpha**2)**s #laplace系数的被积函数
    func=lambda a,b:(b-a)/6*(f1(a)+4*f1((a+b)/2)+f1(b))
    def I(x,y,tol0):
        if abs(func(x,y)-func(x,(x+y)/2)-func((x+y)/2,y))<tol0:
            return func(x,y)
        else:
            return I(x,(x+y)/2,tol0/2)+I((x+y)/2,y,tol0/2)  #若大于误差，则减小步长，递归
    return I(0,2*math.pi,err)/math.pi

if __name__=='__main__':
    
    alpha=np.arange(0.01,0.91,0.01)
    y=[]
    for n in range(0,len(alpha),1):
        y.append(laplace(1,1.5,alpha[n],err))
        print(n)
        print(y)
    
    plt.plot(alpha,y)
    
