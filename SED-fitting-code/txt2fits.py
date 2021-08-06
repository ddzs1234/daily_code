from astropy.io import fits
import numpy as np
#read txt
f=open('/media/ashley/zsash/linux/junior/sed/catalog/TENIS_v1.0.cat','r')
p=f.readlines()
p0=[]
for i in range(0,len(p),1):
	p0.append(p[i])
#column
#objid
#ra
#dec
#flux_J
#flux_Ks
objid=[i.strip().split('$')[1] for i in p0]
ra=[i.strip().split('$')[2] for i in p0]
dec=[i.strip().split('$')[3] for i in p0]
flux_J=[i.strip().split('$')[4] for i in p0]
flux_eJ=[i.strip().split('$')[5] for i in p0]
flux_Ks=[i.strip().split('$')[6] for i in p0]
flux_eKs=[i.strip().split('$')[7] for i in p0]
#print 'objid:',objid
#print type(objid[1]),type(ra[1]),type(dec[1]),type(flux_J[1]),type(flux_Ks[1])
#create a new fits 
col1=fits.Column(name='id',format='I',array=map(int,objid))
col2=fits.Column(name='ra',format='D',array=map(float,ra))
col3=fits.Column(name='dec',format='D',array=map(float,dec))
col4=fits.Column(name='flux_J',format='D',array=map(float,flux_J))
col5=fits.Column(name='flux_eJ',format='D',array=map(float,flux_eJ))
col6=fits.Column(name='flux_Ks',format='D',array=map(float,flux_Ks))
col7=fits.Column(name='flux_eKs',format='D',array=map(float,flux_eKs))
#create a ColDefs(column-definitions) object for al columns:
cols=fits.ColDefs([col1,col2,col3,col4,col5,col6,col7])
#create a new binary table HDU object by using BinTableHDU.from_columns functionthis function returns a BinTableHDU
tbhdu=fits.BinTableHDU.from_columns(cols)
#hdu=fits.PrimaryHDU(n)
prihdr=fits.Header()
prihdr['COMMENT']="Here's some commentary about this FITS file."
prihdu=fits.PrimaryHDU(header=prihdr)
thdulist=fits.HDUList([prihdu,tbhdu])
thdulist.writeto('table.fits')
f.close()
