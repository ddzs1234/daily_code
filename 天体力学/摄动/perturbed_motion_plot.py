# -*- coding: utf-8 -*-
"""
Created on Fri May 18 17:27:36 2018

@author: surface
"""

#绘制fortran程序计算的结果，包括轨道位置，速度以及用长期摄动和普通牛顿运动方程计算的轨道根数变化
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif']=['SimHei']
mpl.rcParams['axes.unicode_minus']=False
position=[]
time1=[]
lab=['x','y','z','x','y','z']
with open("position.txt") as f:
    b=f.readline()
    while b:
        b=b.split(' ')
        b=[float(x) for x in b if x!='']
        if b[0]>=0:
            time1.append(b[0])
            position.append(b[1:])
        else:
            time1.insert(0,b[0])
            position.insert(0,b[1:])
        b=f.readline()
position=np.array(position)
for i in range(6):
    plt.plot(time1,position[:,i])
    plt.ylabel('$'+lab[i]+' \ (m)$')
    plt.xlabel('time (years)')
    plt.savefig('position_'+str(i)+'.png',format='png')
    plt.show()
del position
velocity=[]
with open("velocity.txt") as f:
    b=f.readline()
    while b:
        b=b.split(' ')
        b=[float(x) for x in b if x!='']
        if b[0]>=0:
            velocity.append(b[1:])
        else:
            velocity.insert(0,b[1:])
        b=f.readline()
velocity=np.array(velocity)
for i in range(6):
    plt.plot(time1,velocity[:,i])
    plt.ylabel('$V_{'+lab[i]+'} \ (m \ s^{-1})$')
    plt.xlabel('time (years)')
    plt.savefig('velocity_'+str(i)+'.png',format='png')
    plt.show()
del velocity
elements1=[]
with open("elements1.txt") as f:
    b=f.readline()
    while b:
        b=b.split(' ')
        b=[float(x) for x in b if x!='']
        if b[0]>=0:
            elements1.append(b[1:])
        else:
            elements1.insert(0,b[1:])
        b=f.readline()
elements1=np.array(elements1)
time=[]
elements=[]
with open("elements.txt") as f:
    b=f.readline()
    while b:
        b=b.split(' ')
        b=[float(x) for x in b if x!='']
        if b[0]>=0:
            time.append(b[0])
            elements.append(b[1:])
        else:
            time.insert(0,b[0])
            elements.insert(0,b[1:])
        b=f.readline()
elements=np.array(elements)
plt.plot(time,elements[:,0],'b',label=u'木星长期线性摄动解')
plt.plot(time1,elements1[:,0],'r',linestyle='-.',label=u'木星一般三体问题解')
plt.plot(time,elements[:,4],'k',label=u'土星长期线性摄动解')
plt.plot(time1,elements1[:,2],'y',linestyle='-.',label=u'土星一般三体问题解')
plt.legend(loc="upper right")
plt.ylabel('e')
plt.xlabel('time (years)')
plt.savefig('eccencity'+'.png',format='png')
plt.show()
plt.plot(time,elements[:,1],'b',label=u'木星长期线性摄动解')
plt.plot(time,elements[:,5],'k',label=u'土星长期线性摄动解')
plt.legend(loc="upper right")
plt.ylabel('$pomega \ (°)$')
plt.xlabel('time (years)')
plt.savefig('pomega'+'.png',format='png')
plt.show()
plt.plot(time,elements[:,2],'b',label=u'木星长期线性摄动解')
plt.plot(time1,elements1[:,1],'r',linestyle='-.',label=u'木星一般三体问题解')
plt.plot(time,elements[:,6],'k',label=u'土星长期线性摄动解')
plt.plot(time1,elements1[:,3],'y',linestyle='-.',label=u'土星一般三体问题解')
plt.legend(loc="upper right")
plt.ylabel('$i \ (°)$')
plt.xlabel('time (years)')
plt.savefig('inclination'+'.png',format='png')
plt.show()
plt.plot(time,elements[:,3],'b',label=u'木星长期线性摄动解')
plt.plot(time,elements[:,7],'k',label=u'土星长期线性摄动解')
plt.legend(loc="upper right")
plt.ylabel('$omega \ (°)$')
plt.xlabel('time (years)')
plt.savefig('omega.png',format='png')
plt.show()












