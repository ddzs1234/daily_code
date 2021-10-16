# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 18:35:45 2018

@author: surface
@author: ashley
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


#parameters
err=1e-5
mu=0.001
tend=20000

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
        r1_3=((x0[0,0]+mu)**2+x0[1,0]**2+x0[2,0]**2)**1.5
        r2_3=((x0[0,0]+mu-1)**2+x0[1,0]**2+x0[2,0]**2)**1.5
        omega_x=x0[0,0]-(1-mu)*(x0[0,0]+mu)/r1_3-mu*(x0[0,0]+mu-1)/r2_3
        omega_y=x0[1,0]*(1-(1-mu)/r1_3-mu/r2_3)
        omega_z=-x0[2,0]*((1-mu)/r1_3+mu/r2_3)
        y.extend([omega_x+2*x0[4,0],omega_y-2*x0[3,0],omega_z])
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

'''
def point(x,cj,mu):
    r1=((x+mu)**2)**0.5
    r2=((x+mu-1)**2)**0.5
    return -cj-x**2+(2*(1-mu)/r1)+(2*mu/r2)
'''    


p_x=[]
p_xdot=[]
r_all=[[]]
v_all=[[]]
if __name__=='__main__':
    
    
    t=0

    x0=np.mat([0.499,(3)**0.5/2,0,0,0,0]).T
    Time=[t]  #time
    r=[[x0[0,0],x0[1,0],x0[2,0]]]
    v=[[x0[3,0],x0[4,0],x0[5,0]]] #velociy
    r1=((x0[0,0]+mu)**2+x0[1,0]**2+x0[2,0]**2)**0.5
    r2=((x0[0,0]+mu-1)**2+x0[1,0]**2+x0[2,0]**2)**0.5
    cj0=-2*((x0[3,0]**2+x0[4,0]**2+x0[5,0]**2)/2-(x0[0,0]**2+x0[1,0]**2)/2-(1-mu)/r1-mu/r2)
    CJ=[cj0]
    while t<tend:
        h=0.01
        h,x0_new=RKF78(h,t,x0,err,6)
        x0=x0_new
        t=t+h
        r.append([x0_new[0,0],x0_new[1,0],x0_new[2,0]])
        v.append([x0_new[3,0],x0_new[4,0],x0_new[5,0]])
        Time.append(t)
        r1=((x0_new[0,0]+mu)**2+x0_new[1,0]**2+x0_new[2,0]**2)**0.5
        r2=((x0_new[0,0]+mu-1)**2+x0_new[1,0]**2+x0_new[2,0]**2)**0.5
        cj0=-((x0_new[3,0]**2+x0_new[4,0]**2+x0_new[5,0]**2)-(x0_new[0,0]**2+x0_new[1,0]**2)-2*(1-mu)/r1-2*mu/r2)
        CJ.append(cj0)
    r=np.array(r)
    v=np.array(v)
            
#    ry=np.ndarray.tolist(r[:,1])
#    vy=np.ndarray.tolist(v[:,1])
#        
#        #vx=np.nda
#    test_y=[]
#    for i in range(0,len(ry)):
#        if np.abs(ry[i]-(3)*0.05/2)<0.001 and vy[i]>0 :
#            test_y.append(ry[i])
#            
#        
#    for n in range(0,len(test_y),1):
#        index=ry.index(test_y[n])
#        p_x.append(r[index,0])
#        p_xdot.append(v[index,0])
      
          

    #plt.scatter(r[:,0],r[:,1],s=4)
        #plt.savefig('./L2/L2_orbit.png',format='png')
    #plt.figure()
    CJ=[round(n,4) for n in CJ]
    plt.plot(Time,CJ)
        #plt.savefig('./L2/L2_cj.png')
    
    
        
        #plt.savefig('./L2/L2_p.png',format='png')
























