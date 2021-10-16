#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 14:52:15 2018

@author: ashley
"""

import glob
from astropy.io import fits
import numpy as np
from astropy.table import Table,join,vstack,hstack
from PyAstronomy.pyasl import ftest 


def collect(filename,parameter,fmt): 
    f=fits.open(filename)
    data=f[1].data
    column=[]
    for i in range(0,len(parameter),1):
        c=fits.Column(name=parameter[i],array=np.array(data.field(parameter[i])),format=fmt[i])
        column.append(c)
        
    column=fits.ColDefs(column)
    return column
def collect1(filename,parameter):
    f=fits.open(filename)
    data=f[1].data
    column=[]
    for i in range(0,len(parameter),1):
        c=data.field(parameter[i])
        column.append(c[0])
    return column
    

def add(column_old,column_new,column):
    column_old=fits.BinTableHDU.from_columns(column_old+column_new+column)
    return column_old

if __name__ =='__main__':   
    
#    f1=fits.open('/home/ashley/NewDisk/MaGNA/code/tensorflow_image_classifier-master/src/combined_map_v1.fits',format='fits')
#    f2=fits.open('/home/ashley/NewDisk/MaGNA/code/tensorflow_image_classifier-master/src/combined_flux_v1.fits',format='fits')    
#    
#    data1=Table(f1[1].data)
#    data2=Table(f2[1].data)
#    
#    new=join(data1,data2,keys='PLATEIFU')
#    #map+flux
#    new.write('/home/ashley/NewDisk/MaGNA/code/tensorflow_image_classifier-master/src/compare_v1.fits')

    
    f=fits.open('/home/ashley/NewDisk/MaGNA/code/tensorflow_image_classifier-master/src/compare_v2.fits')
    data=Table(f[1].data)
#    print data
#    print data
    
    data.remove_rows(0)
#    print data
    data.add_row(['10001-12701','0.11149','0.88851','intermediate','0.99303','0.00697','clumpy'])
    data.add_row(['10001-12702','0.99484','0.00516', 'clumpy','1.00000','0.00000','clumpy'])
    data.add_row(['10001-12703','0.86448','0.13552','intermediate','0.99915','0.00085','clumpy'])
    data.add_row(['10001-12704','0.98328','0.01672','clumpy','0.99939','0.00061','clumpy'])
    data.add_row(['10001-12705','0.94145','0.05855','clumpy','0.99924','0.00076','clumpy'])
    
    data.write('/home/ashley/NewDisk/MaGNA/code/tensorflow_image_classifier-master/src/compare_v3.fits')
    
    
    
    
    
    
    
    
    
