#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#Created on Tue Jan 23 08:16:14 2018
#@author: ashley
from __future__ import division
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np 
from math import radians,cos,sin,asin,sqrt
import math
from astropy.io import fits
from astropy.table import Table
import fitsio
from fitsio import FITS,FITSHDR


def lamda2v(x):
    return 3e18/x


def obs2emit(v,z):
    return (1+z)*v


def find_all_index(arr,item):
    return[i for i,a in enumerate(arr) if a==item]
    
    
def nummatch(f,v):
    f_da=[i for i in f if i>0]
    f_fu=[i for i in f if i<=0]
    index=[]
    for i in range(0,len(f_da),1):
        index_fu=f.index(f_da[i])
        index.append(index_fu)#dayuling de zuobiao 
    v_final=[]
    if len(f_fu)>0:
        for i in range(0,len(index),1):
            index1=index[i]
            v_final.append(v[index1])
        return v_final
    else:
        return v
    
    
def mag2flux(m):
    return 10**((m+48.6)/(-2.5))


def flux2vL(f,v,d):
    return 4*math.pi*d*d*f*v


def red2dis(z):
    return 3E8*z*3.08568E24/(6.773999999999998e4)


def dis(ra,dec,ra1,dec1):
    ra,dec,ra1,dec1=map(radians,[ra,dec,ra1,dec1])
    dra=ra-ra1
    ddec=dec-dec1
    a=sin(ddec/2)**2+cos(dec)*cos(dec1)*sin(dra/2)**2
    c=2*asin(sqrt(a))
    return c


def flux2flux(flux):
    f=10**((-73.6/2.5)+math.log10(flux))
    return f

    
    
#add herschel data    
zall=fits.open('/home/ashley/Link_sed/catalog/cdfs.v1.6.9.fits')
zher=fits.open('/home/ashley/Link_sed/catalog/herschel.fits')
ra_her=np.ndarray.tolist(zall[1].data.field('ra'))
dec_her=np.ndarray.tolist(zall[1].data.field('dec'))


#catalog read
main=fits.open('/home/ashley/Link_sed/catalog/maincat.fits')
CANDELS=fits.open('/home/ashley/Link_sed/catalog/CANDELS.GOODSS.F160W.v1.fits')
ECDFS=fits.open('/home/ashley/Link_sed/catalog/ecdfs_BVRdet_Subaru_v11.fits')
TENIS=fits.open('/home/ashley/Link_sed/catalog/table.fits')
GALEX=fits.open('/home/ashley/Link_sed/catalog/CDFS_00-xd-mcat.fits')

#ra_dec information
ra_main=np.ndarray.tolist(main[1].data.field('RA'))
dec_main=np.ndarray.tolist(main[1].data.field('DEC'))
cpcat=np.ndarray.tolist(main[1].data.field('CP_CAT'))
ra_candels=np.ndarray.tolist(CANDELS[1].data.field('ra'))
dec_candels=np.ndarray.tolist(CANDELS[1].data.field('dec'))
ra_ecdfs=np.ndarray.tolist(ECDFS[1].data.field('RA'))
dec_ecdfs=np.ndarray.tolist(ECDFS[1].data.field('DEC'))
ra_tenis=np.ndarray.tolist(TENIS[1].data.field('ra'))
dec_tenis=np.ndarray.tolist(TENIS[1].data.field('dec'))
ra_galex=np.ndarray.tolist(GALEX[1].data.field('alpha_j2000'))
dec_galex=np.ndarray.tolist(GALEX[1].data.field('delta_j2000'))

 
#judge
min_dis=0.5*math.pi/(3600*180)#2.42406840554768e-06

flux_info=[]
#delete VLT data
columns=('id','Z','24_WL','24','24_E','100_WL','100','100_E','160_WL','160','160_E',
         'U_CTIO_WL','U_CTIO','U_CTIO_E',
          'U_VIMOS_WL','U_VIMOS','U_VIMOS_E','F435W_WL','F435W','F435W_E',
          'F606W_WL','F606W','F606W_E',
          'F775W_WL','F775W','F775W_E','F814W_WL','F814W','F814W_E',
          'F850LP_WL','F850LP','F850LP_E','F098M_WL','F098M','F098M_E',
          'F105W_WL','F105W','F105W_E','F125W_WL','F125W','F125W_E',
          'F160W_WL','F160W','F160W_E',
          'Ks_ISAAC_WL','Ks_ISAAC','Ks_ISAAC_E','Ks_HAWKI_WL','Ks_HAWKI','Ks_HAWKI_E',
          '3.6_WL','3.6','3.6_E','4.5_WL','4.5','4.5_E',
          '5.8_WL','5.8','5.8_E','8.0_WL','8.0','8.0_E',
          'J_WL','J','J_E',
          'Ks_WL','Ks','Ks_E','NUV_WL','NUV','NUV_E',
          'FUV_WL','FUV','FUV_E',
          'U38_WL','U38','U38_E',
          'IA427_WL','IA427','IA427_E',
          'IA464_WL','IA_464','IA464_E','IA484_WL','IA484','IA484_E',
          'IA505_WL','IA505','IA505_E','IA527_WL','IA527','IA527_E',
          'IA574_WL','IA574','IA574_E','IA624_WL','IA624','IA624_E',
          'IA679_WL','IA679','IA679_E','IA709_WL','IA709','IA709_E',
          'IA738_WL','IA738','IA738_E','IA768_WL','IA768','IA768_E',
          'IA827_WL','IA827','IA827_E') #wl=A flux=jy 
#,
          #
    
wl_candels=[3734,3722,4317,5918,7693,8047,9055,9851,10550,12486,15370,21605,21463,35508,44960,57245,78840]
wl_tenis=[12481,21338]
wl_galex=[2278,1543]
wl_ecdfs=[3706,4253,4631,4843,5059,5256,5760,6227,6778,7070,7356,7676,8243]
wl_her=[240000,979036.1,1539451.3]######## exact effective lambda

n=1
i=45
tip=0
count=0

while i<1009:
    
   
    min_y=[]
    flux_data=[]

    #ra dec information according to cpcat
    if cpcat[i-1]=='CANDELS':
        ra_main_one=main[1].data[i-1][29]
        dec_main_one=main[1].data[i-1][30]
    if cpcat[i-1]=='GEMS   ':
        ra_main_one=main[1].data[i-1][26]
        dec_main_one=main[1].data[i-1][27]
    if cpcat[i-1]=='WFI    ':
        ra_main_one=main[1].data[i-1][20]
        dec_main_one=main[1].data[i-1][21]
    if cpcat[i-1]=='TENIS  ':
        ra_main_one=main[1].data[i-1][32]
        dec_main_one=main[1].data[i-1][33]
    if cpcat[i-1]=='SEDS   ':
        ra_main_one=main[1].data[i-1][35]
        dec_main_one=main[1].data[i-1][36]
    if cpcat[i-1]=='GOODS-S':
        ra_main_one=main[1].data[i-1][23]
        dec_main_one=main[1].data[i-1][24]
    if cpcat[i-1]=='VLA    ':
        ra_main_one=main[1].data[i-1][38]
        dec_main_one=main[1].data[i-1][39]
    if cpcat[i-1]=='...':
        ra_main_one=ra_main[i-1]
        dec_main_one=dec_main[i-1]

    index_main=[]
    index_candels=[]
    
    distance_main_candels=[]
    distance_main_ecdfs=[]
    distance_main_tenis=[]
    distance_main_galex=[]
    distance_main_irac=[]
    distance_main_her=[]
    for j in range(0,len(ra_candels),1):
        ra_candels_one=ra_candels[j]
        dec_candels_one=dec_candels[j]
        distance_main_candels.append(dis(ra_main_one,dec_main_one,ra_candels_one,dec_candels_one))
    for j in range(0,len(ra_ecdfs),1):
        ra_ecdfs_one=ra_ecdfs[j]
        dec_ecdfs_one=dec_ecdfs[j]
        distance_main_ecdfs.append(dis(ra_main_one,dec_main_one,ra_ecdfs_one,dec_ecdfs_one))
    for j in range(0,len(ra_tenis),1):
        ra_tenis_one=ra_tenis[j]
        dec_tenis_one=dec_tenis[j]
        distance_main_tenis.append(dis(ra_main_one,dec_main_one,ra_tenis_one,dec_tenis_one))
    for j in range(0,len(ra_galex),1):
        ra_galex_one=ra_galex[j]
        dec_galex_one=dec_galex[j]
        distance_main_galex.append(dis(ra_main_one,dec_main_one,ra_galex_one,dec_galex_one))
    for j in range(0,len(ra_her),1):
        ra_her_one=ra_her[j]
        dec_her_one=dec_her[j]
        distance_main_her.append(dis(ra_main_one,dec_main_one,ra_her_one,dec_her_one))
        

    dis_candels=np.min(distance_main_candels)
    dis_ecdfs=np.min(distance_main_ecdfs)
    dis_tenis=np.min(distance_main_tenis)
    index_dis_tenis=distance_main_tenis.index(dis_tenis)
    dis_galex=np.min(distance_main_galex)
    dis_her=np.min(distance_main_her)
    
    

    if main[1].data[i-1][50]>0 and dis_her<=min_dis:#and dis_her <= min_dis:#redshift final
        
        flux_data.append(int(i))
        z_final=main[1].data[i-1][50]
        D=red2dis(z_final)
        flux_data.append(z_final)

        if dis_her <= min_dis:
            index_her=distance_main_her.index(dis_her)
            m=0
            #if  zher[1].data[index_her][2]>0 and zher[1].data[index_her][4]>0 and zher[1].data[index_her][6]>0:
             #   count=count+1
             #   print('i',i)
            
            for j in range(2,7,2):
                if zher[1].data[index_her][j]>0:
                    flux_data.append(wl_her[m])
                    flux_data.append(zher[1].data[index_her][j]*10**(-3))
                    flux_data.append(zher[1].data[index_her][j+1]*10**(-3))
                else:
                    flux_data.append(wl_her[m])
                    flux_data.append(-9999)
                    flux_data.append(-9999)
                m=m+1
        else:
            for j in range(0,3,1):
                flux_data.append(wl_her[j])
                flux_data.append(-9999)
                flux_data.append(-9999)
        
        if dis_candels <= min_dis:
            index_candels=distance_main_candels.index(dis_candels)
            m=0
            for j in range(7,40,2):
                
                if CANDELS[1].data[index_candels][j]>0:
                    flux_data.append(wl_candels[m])
                    flux_data.append(CANDELS[1].data[index_candels][j]*(10**(-6)))
                    flux_data.append(CANDELS[1].data[index_candels][j+1]*(10**(-6)))
                else:
                    flux_data.append(wl_candels[m])
                    flux_data.append(-9999)
                    flux_data.append(-9999)
                m=m+1
        else:
            for j in range(0,17,1):
                flux_data.append(wl_candels[j])
                flux_data.append(-9999)
                flux_data.append(-9999)
         
            
        if dis_tenis <= min_dis:
            index_tenis=distance_main_tenis.index(dis_tenis)
            m=0
            for j in range(3,6,2):
                if TENIS[1].data[index_tenis][j]>0:
                    flux_data.append(wl_tenis[m])
                    flux_data.append(TENIS[1].data[index_tenis][j]*(10**(-6)))
                    flux_data.append(TENIS[1].data[index_tenis][j+1]*(10**(-6)))
                else:
                    flux_data.append(wl_tenis[m])
                    flux_data.append(-9999)
                    flux_data.append(-9999)
                m=m+1
        else:
            for j in range(0,2,1):
                flux_data.append(wl_tenis[j])
                flux_data.append(-9999)
                flux_data.append(-9999)
        
        if dis_galex <= min_dis:
            index_galex=distance_main_galex.index(dis_galex)
            m=0
            for j in range(34,37,2):
                if GALEX[1].data[index_galex][j]>0:
                    flux_data.append(wl_galex[m])
                    flux_data.append(GALEX[1].data[index_galex][j]*(10**(-6)))
                    flux_data.append(GALEX[1].data[index_galex][j+1]*(10**(-6)))
                else:
                    flux_data.append(wl_galex[m])
                    flux_data.append(-9999)
                    flux_data.append(-9999)
                m=m+1
        else:
            for j in range(0,2,1):
                flux_data.append(wl_galex[j])
                flux_data.append(-9999)
                flux_data.append(-9999)

        if dis_ecdfs <= min_dis:
            index_ecdfs=distance_main_ecdfs.index(dis_ecdfs)
            bvr_flux_auto=ECDFS[1].data[index_ecdfs][10]
            f_bvr=ECDFS[1].data[index_ecdfs][12]
            totcor=ECDFS[1].data[index_ecdfs][8]
            if ECDFS[1].data[index_ecdfs][14]>0:
                flux_data.append(wl_ecdfs[0])
                flux_data.append(ECDFS[1].data[index_ecdfs][14]*(10**(-6))*bvr_flux_auto*totcor*1.038*1.185/f_bvr)
                flux_data.append(ECDFS[1].data[index_ecdfs][15]*(10**(-6)))
            else:
                flux_data.append(wl_ecdfs[0])
                flux_data.append(-9999)
                flux_data.append(-9999)
            if ECDFS[1].data[index_ecdfs][34]>0:
                flux_data.append(wl_ecdfs[1])
                flux_data.append(ECDFS[1].data[index_ecdfs][34]*(10**(-6))*bvr_flux_auto*totcor*1.038*1.185/f_bvr)
                flux_data.append(ECDFS[1].data[index_ecdfs][35]*(10**(-6)))
            else:
                flux_data.append(wl_ecdfs[1])
                flux_data.append(-9999)
                flux_data.append(-9999)    
            correct1=[1.032*1.054,1.030*0.934,1.028*0.874,1.027*0.990]    
            m=0
            for j in range(38,45,2):
                
                if ECDFS[1].data[index_ecdfs][j]>0:
                    flux_data.append(wl_ecdfs[m+2])
                    flux_data.append(ECDFS[1].data[index_ecdfs][j]*(10**(-6))*correct1[m]*bvr_flux_auto*totcor/f_bvr)
                    flux_data.append(ECDFS[1].data[index_ecdfs][j+1]*(10**(-6)))
                else:
                    flux_data.append(wl_ecdfs[m+2])
                    flux_data.append(-9999)
                    flux_data.append(-9999)
                m=m+1    
            if ECDFS[1].data[index_ecdfs][48]>0:
                flux_data.append(wl_ecdfs[6])
                flux_data.append(ECDFS[1].data[index_ecdfs][48]*(10**(-6))*bvr_flux_auto*totcor*1.024*0.952/f_bvr)
                flux_data.append(ECDFS[1].data[index_ecdfs][49]*(10**(-6)))
            else:
                flux_data.append(wl_ecdfs[6])
                flux_data.append(-9999)
                flux_data.append(-9999)  
            if ECDFS[1].data[index_ecdfs][52]>0:
                flux_data.append(wl_ecdfs[7])
                flux_data.append(ECDFS[1].data[index_ecdfs][52]*(10**(-6))*bvr_flux_auto*totcor*1.021*0.898/f_bvr)
                flux_data.append(ECDFS[1].data[index_ecdfs][53]*(10**(-6)))
            else:
                flux_data.append(wl_ecdfs[7]/(10**4))
                flux_data.append(-9999)
                flux_data.append(-9999) 
            
            
    
            correct2=[1.020*1.033,1.019*1.001,1.018*0.952,1.017*0.924]
            m=0
            for j in range(56,63,2):
                
                if ECDFS[1].data[index_ecdfs][j]>0:
                    flux_data.append(wl_ecdfs[m+8])
                    flux_data.append(ECDFS[1].data[index_ecdfs][j]*(10**(-6))*bvr_flux_auto*totcor*correct2[m]/f_bvr)
                    flux_data.append(ECDFS[1].data[index_ecdfs][j+1]*(10**(-6)))
                else:
                    flux_data.append(wl_ecdfs[m+8])
                    flux_data.append(-9999)
                    flux_data.append(-9999)
                m=m+1
                
            if ECDFS[1].data[index_ecdfs][66]>0:
                flux_data.append(wl_ecdfs[12])
                flux_data.append(ECDFS[1].data[index_ecdfs][66]*(10**(-6))*bvr_flux_auto*totcor*1.015*1.103/f_bvr)
                flux_data.append(ECDFS[1].data[index_ecdfs][66]*(10**(-6)))
            else:
                flux_data.append(wl_ecdfs[12])
                flux_data.append(-9999)
                flux_data.append(-9999)   
        else:
            for j in range(0,13,1):
                flux_data.append(wl_ecdfs[j])
                flux_data.append(-9999)
                flux_data.append(-9999)
                

        if tip==0:
            update=[flux_data]
        else:
            update.append(flux_data)
            
        #print(i,len(flux_data))

        tip=tip+1
        i=i+1
    else:
        i=i+1
    
update=np.array(update)
t=Table(update,names=columns)
t.write('/home/ashley/Link_sed/AGNfitter-master/data/all_451009.txt',format='ascii')

t.write('/home/ashley/Link_sed/AGNfitter-master/data/all_451009.fits',format='fits')
main.close()
CANDELS.close()
ECDFS.close()
TENIS.close()
GALEX.close()
zall.close()
zher.close()
#print('len(count)',count)  






