# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 17:07:59 2018

@author: surface 
@author: ashley

"""
from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path']='C:/360安全浏览器下载/ffmpeg-20180402-02ae52d-win64-static/bin/ffmpeg.exe'
#import pandas as pd
#import sys

#parameters
G=6.67428e-11
m1=1.989e30*1e8
m2=5.965e24
mu=G*(m1+m2)
tend=3*365*86400#10000*365*86400    #s
AU=149597870700      #meters
err=1e-12*AU
a=1.496e11          #meters
e=0.0167086
clen=25e8
peri=1.471e11     #meters(perihelion)
alpha=23.43*math.pi/180         #Obliquity of the ecliptic(degree)
cos_alpha=math.cos(alpha)
sin_alpha=math.sin(alpha)


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



def generate_coef(t,x,h):
    def diff_func(t,x0):
        y=[]
        y.extend([x0[3,0],x0[4,0],x0[5,0]])
        temp=-mu/(x0[0,0]**2+x0[1,0]**2+x0[2,0]**2)**1.5 
        y.extend([temp*x0[0,0],temp*x0[1,0],temp*x0[2,0]])
        y1=np.mat(y)
        return y1.T
    
    k=np.mat(np.zeros((6,13)))
    for i in range(13):
        k[:,i]=diff_func(t+a0[0,i]*h,x+h*k*b[i].T)
    return k


def RKF78(h,t,x,err,n):
    while True:
        k=generate_coef(t,x,h)
        x0=x+h*k*c.T
        x1=x+h*k*c1.T
        temp=x1-x0
        temp_delta=temp[:3,0].T*temp[:3,0]
        delta=temp_delta[0,0]**0.5/h
        if delta<err:
            break
        else:
            h=h/2
    return (h,x0)


def IsEConservation(x0,x0_new):
    K0=0.5*(x0[3,0]**2+x0[4,0]**2+x0[5,0]**2)-mu/(x0[0,0]**2+x0[1,0]**2+x0[2,0]**2)**0.5
    flag=False
    K0_new=(0.5*(x0_new[3,0]**2+x0_new [4,0]**2+x0_new [5,0]**2)\
             -mu/(x0_new [0,0]**2+x0_new [1,0]**2+x0_new [2,0]**2)**0.5)
    K0_ratio=(K0_new-K0)/K0
    if abs(K0_ratio) <1e-10:
        flag=True
    return flag


def rotate_axes(x0):
    x0_axes=np.mat([x0[0,0],x0[1,0],x0[2,0]])
    rotate_factor=np.mat([[1,0,0],[0,cos_alpha,-sin_alpha],
                         [0,sin_alpha,cos_alpha]])
    return rotate_factor*x0_axes.T


def init():
    line.set_data(r[:,0],r[:,1])
    line.set_3d_properties(r[:,2])
    return line,


def update_lines(num,datalines,lines):
    for line,data in zip(lines,datalines):
        line.set_data(data[0:2,num])
        line.set_3d_properties(data[2,num])
        line.set_marker('o')
        ax.set_title('t='+('%-10.4f' %time0[num])+'yrs')
    return lines


if __name__=='__main__':
    t=0
    x0=np.mat([0,a-clen,0,-30300,0,0]).T
    x0_r=rotate_axes(x0)
    r=[[x0_r[0,0],x0_r[1,0],x0_r[2,0]]]
    v=[[x0[3,0],x0[4,0],x0[5,0]]] #velociy
    Time=[t]                                    #time 
    x0=np.mat([r[0][0],r[0][1],r[0][2],-30300,0,0]).T
    interval1=100
    count=0
    v_all=[np.sqrt(x0[3,0]**2+x0[4,0]**2+x0[5,0]**2)]
    while t<tend:
        count+=1
        h=1500
        h,x0_new=RKF78(h,t,x0,err,6)
        if IsEConservation(x0,x0_new):
            x0=x0_new
            t=t+h
            if count%interval1==0:
                r.append([x0_new[0,0],x0_new[1,0],x0_new[2,0]])
                v.append([x0_new[3,0],x0_new[4,0],x0_new[5,0]])
                v_all.append(np.sqrt(x0_new[3,0]**2+x0_new[4,0]**2+x0_new[5,0]**2))
                Time.append(t)
        else:
            print("WARNING:energy conservation false,please retry")
            break
    fig=plt.figure()
    ax = p3.Axes3D(fig)
    line,=ax.plot([],[],[])
    ax.grid()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim3d([-1.5e11, 1.5e11])
    ax.set_ylim3d([-1.5e11, 1.5e11])
    ax.set_zlim3d([-7e10, 7e10])
    
    r=np.array(r)
    v=np.array(v)
    time0=np.array(Time)/(365*86400)
    '''
    with open('./data.txt','w') as f:
        for i in range(len(r)):
            f.write(('%-10.4f' %time0[i])+'   ')
            for j in range(3):
                f.write(('%-18.10f' %r[i,j])+'   ')
            for j in range(3):
                f.write(('%-18.10f' %v[i,j])+'   ')
            f.write('\n')
    '''
    data=[np.array([r[:,0],r[:,1],r[:,2]])]
    lines=[ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
    ani = animation.FuncAnimation(fig,update_lines,frames=211,fargs=(data,lines),interval=0.1,repeat=True,init_func=init)
    mywriter = animation.FFMpegWriter()
    ani.save('./test_BH.mp4',writer=mywriter)
    plt.show()
    

    
   


































