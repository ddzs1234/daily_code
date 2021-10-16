# -*- coding: utf-8 -*-
"""
Created on Fri May 11 08:10:31 2018

@author: surface
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import math
import time


#parameters
tend=2*365*86400
m0=1.989e30
num=6
g_const=6.67428e-11
AU=149597870700      #meters
err=1e-10
err1=AU*1e-6

eccen=np.array([0.00677323,0.01671022,0.04839266,0.05415060,0.04716771,0.00858587])
inclination=np.array([3.39471,0.00005,1.30530,2.48446,0.76986,1.76917])*np.pi/180
omega=np.array([76.68069,360-11.26064,100.55615,113.71504,74.22988,131.72169])*np.pi/180
pomega=np.array([131.53298,102.94719,14.75385,92.43194,170.96424,44.97135])*np.pi/180
lamda=np.array([181.97973,100.46435,34.40438,49.94432,313.23218,304.88003])*np.pi/180
semi_axis=np.array([0.72333199,1.00000011,5.20336301,9.53707032,19.19126393,30.06896348])*AU
#n1=np.array([210664136.06,129597740.63,10925078.35,4401052.95,1542547.79,786449.21])
#n0=2*math.pi*n1/(360000*360*86400*365)
mass=np.array([48.685e23,59.736e23,18986e23,5684.6e23,868e23,1024.3e23])
n0=(g_const*(m0+mass)/semi_axis**3)**0.5
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
i3=np.mat([[1,0,0],[0,1,0],[0,0,1]])
#semi_axis=[5.20336301,9.53707032]
alpha=np.zeros((num,num))
alpha=np.mat(alpha)
for i in range(num):
    for j in range(i+1,num):
        alpha[i,j]=semi_axis[i]/semi_axis[j]
alpha=alpha+alpha.T
c00=1.25*np.array(alpha)*(1-0.125*np.array(alpha)**2)


def b_alpha(x):
    return 3*x*(1+15/8*x**2+175/64*x**4+3675/1024*x**6+72765/16384*x**8)

c_i=np.zeros((num,num))
for i in range(num):
    for j in range(num):
        c_i[i,j]=0.25*n0[i]*mass[j]/(m0+mass[i])*alpha[i,j]**2*b_alpha(alpha[i,j])

A1=np.zeros((num,num))
for i in range(num):
    A1[i,i]=sum(c_i[i,:])
A2=-c00*c_i
A=A1+A2
B=A1-c_i
tran_mat=np.zeros((4*num,4*num))
tran_mat[0:num,num:2*num]=A
tran_mat[num:2*num,0:num]=-A
tran_mat[2*num:3*num,3*num:4*num]=B
tran_mat[3*num:4*num,2*num:3*num]=-B
tran_mat=np.mat(tran_mat)

def generate_coef(t,x,h,n):
    def diff_func(t,x0):
        return tran_mat*x0
    
    k=np.mat(np.zeros((n,13)))
    for i in range(13):
        k[:,i]=diff_func(t+a0[0,i]*h,x+h*k*b[i].T)
    return k


def RKF78(h,t,x,err,n):
        k=generate_coef(t,x,h,n)
        x0=x+h*k*c.T
        x1=x+h*k*c1.T
        temp=x1-x0
        temp_delta=temp.T*temp
        delta=temp_delta[0,0]**0.5/h
        if delta>err:
            x00=RKF78(h/2,t,x,err,n)
            t=t+h/2
            x0=RKF78(h/2,t,x00,err,n)
        return x0

def generate_coef_1(t,x,h,n):
    def diff_func_1(t,x0):
        mat0=np.mat(np.zeros((6*num,6*num)))
        for i in range(num):
            mat0[3*i:3*i+3,3*num+3*i:3*num+3*i+3]=i3
        r_3=[]
        for i in range(num):
            r_3.append((x0[3*i:3*i+3,0].T*x0[3*i:3*i+3,0])[0,0]**1.5)
        rij_3=np.zeros((num,num))
        for i in range(num):
            for j in range(i):
                rij_3[i,j]=((x0[3*i:3*i+3,0]-x0[3*j:3*j+3,0]).T*(x0[3*i:3*i+3,0]-x0[3*j:3*j+3,0]))[0,0]**1.5
        rij_3=rij_3+rij_3.T
        for i in range(num):
            for j in range(num):
                if i==j:
                    mat0[3*num+3*i:3*num+3*i+3,3*i:3*i+3]=(-g_const*(m0+mass[i])/r_3[i]-g_const*\
                        sum([mass[k]/rij_3[i,k] for k in range(num) if k!=i]))*i3
                else:
                    mat0[3*num+3*i:3*num+3*i+3,3*j:3*j+3]=g_const*mass[j]*(1/rij_3[i,j]-1/r_3[j])*i3
        return mat0*x0
    
    k=np.mat(np.zeros((6*num,13)))
    for i in range(13):
        k[:,i]=diff_func_1(t+a0[0,i]*h,x+h*k*b[i].T)
    return k


def RKF78_1(h,t,x,err,n):
        k=generate_coef_1(t,x,h,n)
        x0=x+h*k*c.T
        x1=x+h*k*c1.T
        temp=x1-x0
        temp_delta=temp[:3*num,0].T*temp[:3*num,0]
        delta=temp_delta[0,0]**0.5/h
        if delta>err:
            x00=RKF78_1(h/2,t,x,err,n)
            t=t+h/2
            x0=RKF78_1(h/2,t,x00,err,n)
        return x0


def M2E(e,m):
    sigma=1e-300
    delta=1
    e0=0
    e1=0
    while delta>sigma:
        e1=e0-(e0-e*math.sin(e0)-m)/(1-e*math.cos(e0))
        delta=abs(e1-e0)
        e0=e1
    return e1


def rx(s):
    return np.mat([[1,0,0],[0,math.cos(s),math.sin(s)],[0,-math.sin(s),math.cos(s)]])


def ry(s):
    return np.mat([[math.cos(s),0,-math.sin(s)],[0,1,0],[math.sin(s),0,math.cos(s)]])

def rz(s):
    return np.mat([[math.cos(s),math.sin(s),0],[-math.sin(s),math.cos(s),0],[0,0,1]])


if __name__=='__main__':
    start=time.time()
    h0=eccen*np.sin(pomega)
    k0=eccen*np.cos(pomega)
    p0=inclination*np.sin(omega)
    q0=inclination*np.cos(omega)
    vector0=np.zeros((1,num*4))
    vector0[0,0:num]=h0
    vector0[0,num:2*num]=k0
    vector0[0,2*num:3*num]=p0
    vector0[0,3*num:4*num]=q0
    vector0=np.mat(vector0)
    vector0=vector0.T
    vector1=[]
    target=1
    vector1.append([vector0[target,0],vector0[target+num,0],\
                    vector0[target+2*num,0],vector0[target+3*num,0]])
    t=0
    h=86400
    Time=[t]
    while t<tend:
        vector0_new=RKF78(h,t,vector0,err,4*num)
        vector0=vector0_new
        t=t+h
        vector1.append([vector0[target,0],vector0[target+num,0],\
                        vector0[target+2*num,0],vector0[target+3*num,0]])
        Time.append(t)
    vector1=np.array(vector1)
    elements=np.zeros((len(Time),4))
    elements[:,0]=(vector1[:,0]**2+vector1[:,1]**2)**0.5
    temp=np.arccos(vector1[:,1]/elements[:,0])
    elements[:,1]=(temp+(np.pi-temp)*np.sign(vector1[:,0])*(np.sign(vector1[:,0])-1))*180/np.pi
    elements[:,2]=(vector1[:,2]**2+vector1[:,3]**2)**0.5
    temp=np.arccos(vector1[:,3]/elements[:,2])
    elements[:,3]=(temp+(np.pi-temp)*np.sign(vector1[:,2])*(np.sign(vector1[:,2])-1))*180/np.pi
    elements[:,2]=elements[:,2]*180/np.pi
    Time=np.array(Time)/86400/365
    end=time.time()
    print("calculation 1 end")
    print("time used: "+str(end-start))
    #下面是用一般的牛顿运动方程求解
    start=time.time()
    M2E0=np.frompyfunc(M2E,2,1)
    m_temp=[]
    for m in lamda-pomega:
        if m<0:
            m=m+2*np.pi
        m_temp.append(m)
    e0=M2E0(eccen,m_temp)
    x0=np.zeros((6*num,1))
    for i in range(num):  #前一半为位置坐标
        temp=np.mat([[semi_axis[i]*(math.cos(e0[i])-eccen[i])],\
                      [semi_axis[i]*(1-eccen[i]**2)**0.5*math.sin(e0[i])],[0]])
        temp1=np.array(rz(-omega[i])*rx(-inclination[i])*rz(omega[i]-pomega[i])*temp)
        x0[3*i:3*i+3,0]=temp1[:,0]
    for i in range(num):  #后一半为速度坐标
        r=sum(x0[3*i:3*i+3,0]**2)**0.5
        temp=np.mat([[-(semi_axis[i])**2*n0[i]/r*math.sin(e0[i])],\
                      [(semi_axis[i])**2*n0[i]/r*(1-eccen[i]**2)**0.5*math.cos(e0[i])],[0]])
        temp1=np.array(rz(-omega[i])*rx(-inclination[i])*rz(omega[i]-pomega[i])*temp)
        x0[3*num+3*i:3*num+3*i+3,0]=temp1[:,0]
    x0=np.mat(x0) #得到初始时刻各行星的位置和速度
    t=0
    h=864
    Time1=[t]
    position=[]
    position.append([x0[3*target,0],x0[3*target+1,0],x0[3*target+2,0]])
    velocity=[]
    velocity.append([x0[3*num+3*target,0],x0[3*num+3*target+1,0],x0[3*num+3*target+2,0]])
    while t<tend:
        x0_new=RKF78_1(h,t,x0,err1,6*num)
        x0=x0_new
        t=t+h
        position.append([x0[3*target,0],x0[3*target+1,0],x0[3*target+2,0]])
        velocity.append([x0[3*num+3*target,0],x0[3*num+3*target+1,0],x0[3*num+3*target+2,0]])
        Time1.append(t)
    position=np.array(position)
    fig=plt.figure()
    ax=p3.Axes3D(fig)
    ax.plot(position[:,0],position[:,1],position[:,2])
    plt.show()
    velocity=np.array(velocity)
    elements1=np.zeros((len(Time1),4))
    r0=(position[:,0]**2+position[:,1]**2+position[:,2]**2)**0.5
    v0_2=velocity[:,0]**2+velocity[:,1]**2+velocity[:,2]**2
    semi=1/(2/r0-v0_2/(g_const*(m0+mass[target])))
    n2=(g_const*(m0+mass[target])/semi**3)**0.5
    h=((position[:,1]*velocity[:,2]-position[:,2]*velocity[:,1])**2+\
       (position[:,2]*velocity[:,0]-position[:,0]*velocity[:,2])**2+\
       (position[:,0]*velocity[:,1]-position[:,1]*velocity[:,0])**2)**0.5
    elements1[:,0]=(1-h**2/(g_const*(m0+mass[target])*semi))**0.5
    hz=position[:,0]*velocity[:,1]-position[:,1]*velocity[:,0]
    elements1[:,2]=np.arccos(hz/h)*180/np.pi
    del hz,h,v0_2
    hx=position[:,1]*velocity[:,2]-position[:,2]*velocity[:,1]
    hy=position[:,2]*velocity[:,0]-position[:,0]*velocity[:,2]
    temp=np.arccos(-hy/(hx**2+hy**2)**0.5)
    elements1[:,3]=(temp+(np.pi-temp)*np.sign(hx)*(np.sign(hx)-1))*180/np.pi
    u=(1-r0/semi)/elements1[:,0]
    v=(position[:,0]*velocity[:,0]+position[:,1]*velocity[:,1]+\
       position[:,2]*velocity[:,2])/(n2*semi**2*elements1[:,0])
    temp=np.arccos((v*position[:,2]/(r0*(1-elements1[:,0]**2))+\
                   (u-elements1[:,0])*velocity[:,2]/(n2*semi*\
                   (1-elements1[:,0]**2)**0.5))/np.sin(elements1[:,2]))
    temp1=u*position[:,2]/r0-v*velocity[:,2]/(n2*semi)
    elements1[:,1]=(temp+(np.pi-temp)*np.sign(temp1)*(np.sign(temp1)-1))*180/np.pi
    end=time.time()
    print("calculation 2 end")
    print("time used: "+str(end-start))
    Time1=np.array(Time1)/86400/365
    plt.plot(Time,elements[:,0],'b',Time1,elements1[:,0],'r')
    plt.show()
    #plt.plot(Time,elements[:,1],'b',Time1,elements1[:,1],'r')
    #plt.show()
    plt.plot(Time,elements[:,2],'b',Time1,elements1[:,2],'r')
    plt.show()
    #plt.plot(Time,elements[:,3],'b',Time1,elements1[:,3],'r')
    #plt.show()






















