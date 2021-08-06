#!/usr/bin/python
import numpy as np
import math
from astropy.io import fits
from collections import Counter
import matplotlib.pyplot as plt
hdu=fits.open('/home/ashley/xray/sdss_zs.fit')
f=open('/home/ashley/xray/cxo-sdssmatch.file','r')
p=f.readlines()
p0=[]
print len(p)
for i in range(3,len(p)-1,1):
	p0.append(p[i])
objid=[n.strip().split('$')[3] for n in p0]
flux=[n.strip().split('$')[5] for n in p0]
print 'objid:',objid
print  'fits:',hdu[1].data.field(0)
s=map(str,hdu[1].data.field(0))
print 's:',s
#result=[i for i in objid if i in hdu[1].data.field(7)]
result=list(set(objid).intersection(set(map(str,hdu[1].data.field(0)))))
print 'result:',result
print 'count:',len(result)
c_sdss=[]
c_csc=[]
for i in range(0,len(result),1):
	c_sdss.append(s.index(result[i]))
	c_csc.append(objid.index(result[i]))
print 'c_sdss:',c_sdss
print 'c_csc:',c_csc
z=[]
flux2=[]
for i in range(0,len(c_sdss),1):
	z.append(hdu[1].data[c_sdss[i]][7])
	flux2.append(flux[c_csc[i]])
print 'z:',z

v=[i*3E8 for i in z]
D=[i*308568E16/67.73999999999998 for i in v]#cm
print 'Distance:',D
print 'flux2:',flux2
L=[]
for i in range(0,len(D),1):
	L.append(float(flux2[i])*4*math.pi*D[i]*D[i])
print 'L:',L

plt.scatter(D,L)
plt.xlabel('Distance'+' '+'(cm)')
plt.ylabel('L'+' '+'(erg'+' '+'s^-1)')
plt.yscale('log')
plt.xscale('log')
#plt.title('hushu')
plt.savefig('cxc_sdss',format='png')



plt.close()
hdu.close()
f.close()
