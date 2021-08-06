from astropy.io import fits
hdu=fits.open('/home/ashley/sed/maincat.fits')#read
#from astropy.cosmology import Planck15 as cosmo
#cosmo.H(0) 67.73999999999998 km / (Mpc s) 1MPC=308568E16cm
import numpy as np
import math
import matplotlib.pyplot as plt
f_eff=[4.56E15,3.52E15,3.33E15,2.4E15,1.37E15,8.3E11,1.4E10]#HZ
i=0#i:f_eff;
j=0#j:id
m=1
c=3E8
f_emitlog=[]
fL=[]#zongzuobiao
colors=['b','g','r','c','m','k','0.75']


#read fig7_R06SED
f=open('/home/ashley/sed/fig7_R06SED.dat','r')
p=f.readlines()
p0=[]
for n in range(1,len(p),1):
    p0.append(p[n])
lognu=[n.strip().split('    ')[0] for n in p0]
lognulnnu=[n.strip().split('    ')[1] for n in p0]

linestyle=['dashed']
while j<2:
    plt.figure(m)
    fL=[]
    f_emitlog=[]
    if hdu[1].data[j][50]==-1 or hdu[1].data[j][50]==0:#redshift null
        j=j+1
    else:
        z_final=hdu[1].data[j][50]#redshift
        
        print(z_final)
        D=c*z_final*308568E16/67.73999999999998#distance:cm
        f_emit=[4.56E15*(1+z_final),3.52E15*(1+z_final),3.33E15*(1+z_final),2.4E15*(1+z_final),1.37E15*(1+z_final),8.3E11*(1+z_final),1.4E10*(1+z_final)]
        if hdu[1].data[j][22]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[0]*(10**((hdu[1].data[j][22]+48.6)/(-2.5)))))#WFI
            f_emitlog.append(math.log10(4.56E15*(1+z_final)))
        if hdu[1].data[j][25]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[1]*(10**((hdu[1].data[j][25]+48.6)/(-2.5)))))#GOODS-S
            f_emitlog.append(math.log10(3.52E15*(1+z_final)))
        if hdu[1].data[j][28]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[2]*(10**((hdu[1].data[j][28]+48.6)/(-2.5)))))#GEMS
            f_emitlog.append(math.log10(3.33E15*(1+z_final)))
	       
        if hdu[1].data[j][31]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[3]*(10**((hdu[1].data[j][31]+48.6)/(-2.5)))))#CANDELS_F125W
            f_emitlog.append(math.log10(2.4E15*(1+z_final)))
	         
        if hdu[1].data[j][34]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[4]*(10**((hdu[1].data[j][34]+48.6)/(-2.5)))))#TENIS
            f_emitlog.append(math.log10(1.37E15*(1+z_final)))
	       
        if hdu[1].data[j][37]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[5]*(10**((hdu[1].data[j][37]+48.6)/(-2.5)))))#SEDS
            f_emitlog.append(math.log10(8.3E11*(1+z_final)))
	         
        if hdu[1].data[j][40]!=-1:
            fL.append(math.log10(4*math.pi*D*D*f_emit[6]*(10**((hdu[1].data[j][40]+48.6)/(-2.5)))))#VLA
            f_emitlog.append(math.log10(1.4E10*(1+z_final)))
            
        print (fL,f_emitlog)
        plt.axis([0.9*np.min(f_emitlog),19.5,0.98*np.min(fL),47])
	    
        plt.scatter(f_emitlog,fL,c=colors,marker='v')
        plt.plot(lognu,lognulnnu,linestyle='dashed')
        plt.xlabel('log'+' '+'$v_{rest}$'+'(Hz)')
        plt.ylabel('log'+' '+'$v$ $L_v$'+'($erg$ $s^{-1}$)')
            
        plt.text(0.93*np.min(13),0.99*np.min(fL)+0.4,'id:'+' '+str(hdu[1].data[j][0]))#id
        plt.text(0.93*np.min(13),0.99*np.min(fL),'z_final:'+' '+str(hdu[1].data[j][50]))#z_final
        plt.text(1.25*np.min(13),0.99*np.min(fL),'VLA',color='0.75')
        plt.text(1.25*np.min(13),0.99*np.min(fL)+0.4,'SEDS',color='k')  #np.min(FL)wrong
        plt.text(1.25*np.min(13),0.99*np.min(fL)+0.8,'TENIS',color='m')
        plt.text(1.25*np.min(13),0.99*np.min(fL)+1.2,'CANDELS',color='c')
        plt.text(1.25*np.min(13),0.99*np.min(fL)+1.6,'GEMS',color='r')
        plt.text(1.25*np.min(13),0.99*np.min(fL)+2.0,'GOODS-S',color='g')
        plt.text(1.25*np.min(13),0.99*np.min(fL)+2.4,'WFI',color='b')
        picname='/home/ashley/sed/sed_maincat/'+str(hdu[1].data[j][0])+'.png'
        plt.savefig(picname,format='png',dpi=100)
        j=j+1
        m=m+1
        plt.close()
        
hdu.close()
