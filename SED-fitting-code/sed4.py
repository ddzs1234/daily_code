# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#!/usr/lib/python
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math
main=fits.open('/home/ashley/sed/maincat.fits')
CANDELS=fits.open('/home/ashley/sed/main_CANDELS.fits')
ECDFS=fits.open('/home/ashley/sed/main_ECDFS.fits')
TENIS=fits.open('/home/ashley/sed/main_TENIS.fits')
#effective frequency
fmain=[4.56E14,3.52E14,3.33E14,2.4E14,1.37E14,8.3E13,1.4E9]#hz
#有效波长 for each catalog
fi=open('/home/ashley/sed/有效波长2014','r')
pi=fi.readlines()
lamda_candels1=[]
lamda_ecdfs1=[]
lamda_tenis1=[]
for i in range(2,len(pi),1):
	if i<=19:
		lamda_candels1.append(pi[i])	
	elif i>19 and i<=51:
		lamda_ecdfs1.append(pi[i])
	elif 52<=i<=53:
		lamda_tenis1.append(pi[i])
lamda_candels=[n.strip().split('\t')[1] for n in lamda_candels1]
lamda_ecdfs=[n.strip().split('\t')[1] for n in lamda_ecdfs1]
lamda_tenis=[n.strip().split('\t')[1] for n in lamda_tenis1]
print(lamda_candels,lamda_ecdfs,lamda_tenis)

#read fig7_R06SED
f=open('/home/ashley/sed/fig7_R06SED.dat','r')
p=f.readlines()
p0=[]
for n in range(1,len(p),1):
    p0.append(p[n])
lognu=[n.strip().split('    ')[0] for n in p0]
lognulnnu=[n.strip().split('    ')[1] for n in p0]
#自定义函数 波长转换为频率 
def lamda2v(x):
    return 3e18/x
fcandels=[lamda2v(i) for i in map(float,lamda_candels)]
fecdfs=[lamda2v(i) for i in map(float,lamda_ecdfs)]
ftenis=[lamda2v(i) for i in map(float,lamda_tenis)]
print('fcandesl',fcandels,fecdfs,ftenis)
#自定义函数 v_rest
def obs2emit(v,z):
    return (1+z)*v
#对于某个AGN在main，分别从candles,ecdfs,tenis match catalog查找

#自定义函数查找某元素在list所有的下标
def find_all_index(arr,item):
    return[i for i,a in enumerate(arr) if a==item]
#自定义函数flux和f——eff的对应：
#hanshuxiecuole
def nummatch(f,v):
    f_da=[i for i in f if i>0]
    #huoqu xiaoyu 0 de zuobiao
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
#自定义函数mag-L
def mag2flux(m):
    return 10**((m+48.6)/(-2.5))
def flux2vL(f,v,d):
    return 4*math.pi*d*d*f*v
def red2dis(z):
    return 3E8*z*308568E16/67.73999999999998
m=1
i=1#id in main

while i<1009:
    plt.figure(m)    
    #获得redshift
    if main[1].data[i-1][50]>0:#redshift final
        z_final=main[1].data[i-1][50]#redshift
        D=red2dis(z_final)
        #main的mag数据
        mag=[]
        for j in range(22,42,3):
            mag.append(main[1].data[i-1][j])
        mag_exist=[n for n in mag if n!= -1]
        flux_main=[mag2flux(n) for n in mag_exist]
        v_main=nummatch(mag,fmain)
        vLmai=[]
        for j in range(0,len(flux_main),1):
            v1mai=v_main[j]
            fluxmai=flux_main[j]
            vLmai.append(flux2vL(fluxmai,v1mai,D))
        logvmai=[math.log10(n) for n in v_main]
        logvLmai=[math.log10(n) for n in vLmai]
        star =plt.scatter(logvmai,logvLmai,marker='*',label='maincat')
        plt.legend(loc=1)
        flux_candels0=[]
        flux_ecdfs0=[]
        flux_tenis0=[] 
        v_all=[]
        flux_all=[]
        #main的id 是否在candels,如果有就获取flux数据：
        #注意flux判断>0：
        #qiyu sange  biaode flux danwei dou buyiyang 
        if i in CANDELS[1].data.field('ID_1'):
            idcandels=np.ndarray.tolist(CANDELS[1].data.field('ID_1'))#id 在candels的第几行
            idindex=idcandels.index(i)
            #flux数据在candels从81列到113列，2递增
            flux_candels=[]
            for j in range(80,114,2):
                flux_candels.append(CANDELS[1].data[idindex][j])#所有的flux
            flux_candels0=[1e-29*n for n in flux_candels if n>0]
            v_candels=nummatch(flux_candels,fcandels)
            if len(flux_candels0)>0:#所有flux>0
                #frequency>0
                v_canrest=[]
                for j in range(0,len(v_candels),1):
                    v1can=obs2emit(v_candels[j],z_final)
                    v_canrest.append(v1can)#rest v>0
                vLcan=[]
                for n in range(0,len(flux_candels0),1):
                    vcan=v_canrest[n]
                    fluxcan=flux_candels0[n]
                    vLcan.append(flux2vL(fluxcan,vcan,D))
                logvcan=[math.log10(n) for n in v_canrest]
                logvLcan=[math.log10(n) for n in vLcan]
                hline =plt.scatter(logvcan,logvLcan,marker='_',label='CANDELS')
                plt.legend(loc=1)
        if i in ECDFS[1].data.field('ID'):
            idecdfs=np.ndarray.tolist(ECDFS[1].data.field('ID'))
            idindex1=idecdfs.index(i)
            flux_ecdfs=[]
            #flux从82到144,递增2
            for j in range(81,144,2):
                flux_ecdfs.append(ECDFS[1].data[idindex1][j])#所有的flux
            flux_ecdfs0=[1e-29*n for n in flux_ecdfs if n>0]#所有flux>0
            v_ecdfs=nummatch(flux_ecdfs,fecdfs)
            if len(flux_ecdfs0)>0:
               
                v_ecdrest=[]
                for j in range(0,len(v_ecdfs),1):
                    v1ecd=obs2emit(v_ecdfs[j],z_final)
                    v_ecdrest.append(v1ecd)
                vLecd=[]
                for n in range(0,len(flux_ecdfs0),1):
                    vecd=v_ecdrest[n]
                    fluxecd=flux_ecdfs0[n]
                    vLecd.append(flux2vL(fluxecd,vecd,D))
                logvecd=[math.log10(n) for n in v_ecdrest]
                logvLecd=[math.log10(n) for n in vLecd]
                vline =plt.scatter(logvecd,logvLecd,marker='|',label='ECDFS')
                plt.legend(loc=1)

        if  i in TENIS[1].data.field('ID_1'):
            idtenis=np.ndarray.tolist(TENIS[1].data.field('ID_1'))
            idindex2=idtenis.index(i)
            flux_tenis=[]
            #flux从82到144,递增2
            for j in range(76,78,1):
                flux_tenis.append(TENIS[1].data[idindex2][j])#所有的flux
            flux_tenis0=[1e-29*n for n in flux_tenis if n>0]#所有flux>0
            if len(flux_tenis0)>0:
                v_tenis=nummatch(flux_tenis,ftenis)#v_tenis>0
                v_tenrest=[]
                for j in range(0,len(flux_tenis0),1):
                    v1ten=obs2emit(v_tenis[j],z_final)
                    v_tenrest.append(v1ten)
                vLten=[]
                for n in range(0,len(flux_tenis0),1):
                    vten=v_tenrest[n]
                    fluxten=flux_tenis0[n]
                    vLten.append(flux2vL(fluxten,vten,D))
                logvten=[math.log10(n) for n in v_tenrest]
                logvLten=[math.log10(n) for n in vLten]
                plus =plt.scatter(logvten,logvLten,marker='+',label='TEINS')
                plt.legend(loc=0)
        plt.plot(lognu,lognulnnu,linestyle='dashed')
        
        plt.text(12.3,46,'id:'+' '+str(main[1].data[i-1][0]))
        plt.text(12.3,45.5,'z_final:'+' '+str(main[1].data[i-1][50]))
        
        plt.xlabel('log'+' '+'$v_{rest}$'+'(Hz)')
        plt.ylabel('log'+' '+'$v$ $L_v$'+'($erg$ $s^{-1}$)')
        picname='/home/ashley/sed/sed4/'+str(main[1].data[i-1][0])+'.png'
        plt.savefig(picname,format='png',dpi=100)
        i=i+1
        m=m+1
        plt.close()
    else:
        i=i+1
    
main.close()
CANDELS.close()
ECDFS.close()
TENIS.close()
       
