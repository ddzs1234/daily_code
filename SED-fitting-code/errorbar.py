#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 13:49:44 2017

@author: ashley
"""

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
'''
import sys
import os
scriptpath='/mnt/wwn-0x5001b444a4ec7c48/xray/code/project/xray/variance.py'
sys.path.append(os.path.abspath(scriptpath))
import variance
'''
file2=fits.open('/mnt/wwn-0x5001b444a4ec7c48/xray/data/dataforyang2.fits')
flux=np.ndarray.tolist(file2[1].data.field('flux_aper_b'))
flux_hilim=np.ndarray.tolist(file2[1].data.field('flux_aper_hilim_b'))
flux_lolim=np.ndarray.tolist(file2[1].data.field('flux_aper_lolim_b'))
name=np.ndarray.tolist(file2[1].data.field('name'))
acis_num1=np.ndarray.tolist(file2[1].data.field('acis_num'))
time=np.ndarray.tolist(file2[1].data.field('gti_mjd_obs'))

i1=46474
#i1=15885
#i2=25691
#i3=46474 CXO J141916.7+525043

acis1=acis_num1[i1]


flux1=[]


flux_hi1=[]


flux_lo1=[]

mjd1=[]


for n in range(i1,acis1+i1,1):
    flux1.append(flux[n])
    flux_hi1.append(flux_hilim[n])
    flux_lo1.append(flux_lolim[n])
    mjd1.append(time[n])
flux1_err1=list(map(lambda x:x[0]-x[1],zip(flux_hi1,flux1)))
flux1_err2=list(map(lambda x:x[0]-x[1],zip(flux1,flux_lo1)))

mjd_detail=[]
mjd_detail1=[]
for n in range(0,len(mjd1),1):
    if 53600<=mjd1[n]<=53700:
        mjd_detail.append(mjd1[n])
    elif mjd1[n]>=54700:
        mjd_detail1.append(mjd1[n])
index_info=[]
index_info1=[]
flux_detail=[]
flux_detail1=[]
flux_detailerr1=[]
flux_detailerr2=[]
flux_detail1err1=[]
flux_detail1err2=[]

for n in range(0,len(mjd_detail),1):
    index=mjd1.index(mjd_detail[n])
    flux_detail.append(flux1[index])
    flux_detailerr1.append(flux1_err1[index])
    flux_detailerr2.append(flux1_err2[index])
    index_info.append(index)

for n in range(0,len(mjd_detail1),1):
    index=mjd1.index(mjd_detail1[n])
    flux_detail1.append(flux1[index])
    flux_detail1err1.append(flux1_err1[index])
    flux_detail1err2.append(flux1_err2[index])
    index_info1.append(index)



#plt.scatter(mjd1,flux1,s=4)
plt.errorbar(mjd1,flux1,yerr=flux1_err1,uplims=True,lolims=False,ms=4,fmt='o')
plt.errorbar(mjd1,flux1,yerr=flux1_err2,lolims=True,uplims=False,ms=4,fmt='o')
plt.title('CXO J141916.7+525043_1')
plt.xlabel('gti-mjd')
plt.ylabel('flux')

sub_axes = plt.axes([.3, .3, .25, .25])
sub_axes.errorbar(mjd_detail,flux_detail,yerr=flux_detailerr1,uplims=True,lolims=False,ms=4,fmt='o')
sub_axes.errorbar(mjd_detail,flux_detail,yerr=flux_detailerr2,lolims=True,uplims=False,ms=4,fmt='o')

sub_axes1 = plt.axes([.5, .5, .25, .25])
sub_axes1.errorbar(mjd_detail1,flux_detail1,yerr=flux_detail1err1,uplims=True,lolims=False,ms=4,fmt='o')
sub_axes1.errorbar(mjd_detail1,flux_detail1,yerr=flux_detail1err2,lolims=True,uplims=False,ms=4,fmt='o')


plt.yscale('log')
#plt.xscale('log')

    
plt.savefig('/mnt/wwn-0x5001b444a4ec7c48/xray/result/CXOJ141916.7+525043_part1.png',format='png')
file2.close()
    
    
    
    
    
    