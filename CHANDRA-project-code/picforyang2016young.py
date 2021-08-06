#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:40:55 2017

@author: ashley
"""

#=======draw the same pic from <yang>

from astropy.io import fits
import matplotlib.pyplot as plt

import sys
import os
scriptpath='/media/ashley/KINGSTON/variance.py'
sys.path.append(os.path.abspath(scriptpath))
import variance

import numpy as np

file2=fits.open('/media/ashley/KINGSTON/dataforyang2.fits')
flux=np.ndarray.tolist(file2[1].data.field('flux_aper_b'))
flux_hilim=np.ndarray.tolist(file2[1].data.field('flux_aper_hilim_b'))
flux_lolim=np.ndarray.tolist(file2[1].data.field('flux_aper_lolim_b'))
name=np.ndarray.tolist(file2[1].data.field('name'))
acis_num1=np.ndarray.tolist(file2[1].data.field('acis_num'))
time=np.ndarray.tolist(file2[1].data.field('gti_mjd_obs'))

max_L=[]
dura=[]
for i in range(0,len(variance.name_2),1):
    index=name.index(variance.name_2[i])
    acis_num=variance.acis_2[i]
    red=variance.redshift_2[i]
    
    dist=red*3E8*308568E16/70
    
    tmp_flux=[]
    tmp_time=[]
    for j in range(index,index + acis_num,1):
        if flux[j]==0.0:
            tmp_flux.append(flux_hilim[j])
            tmp_time.append(time[j])
        elif flux[j]!=0.0 and flux_lolim[j]=='        ':
            tmp_flux.append(flux[j])
            tmp_time.append(time[j])
        else:
            tmp_flux.append(flux[j])
            if flux_hilim[j]-flux[j]>=flux[j]-flux_lolim[j]:
                
                tmp_time.append(time[j])
            else:
                tmp_time.append(time[j])
    max_flux=np.max(tmp_flux)
    max_time=np.max(tmp_time)
    min_time=np.min(tmp_time)
    duratime=(max_time-min_time)*3600
    dura.append(duratime)
    max_l=max_flux*dist**2*4*np.pi
    max_L.append(max_l)

logmax_L=[np.log10(n) for n in max_L]
log_time=[np.log10(n) for n in dura]
plt.figure()
plt.scatter(log_time,logmax_L,s=6,facecolor='none',edgecolor='g')
plt.xlabel('Max log'+'$\Delta$'+'$t$'+'(s)')
plt.ylabel('Max log'+'$L_{x}$'+'(0.5-7.0 kev)'+'($erg$ $s^{-1}$)')
plt.savefig('/media/ashley/KINGSTON/pic1_yang.pdf',dpi=200,format='pdf')

file2.close()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
