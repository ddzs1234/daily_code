"""
unit convert: 
1. WISE: mag 2 flux in Jy
2. VLA: other radio tele: mom0 2 SB_atomgas
"""
import numpy as np

def WISE2flux(magzp,phot):
    """
    https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html
    flux in Jy
    """
    mag    = -2.5*np.log10(phot)+mag_zp
    flux   = 10**(-mag/2.5)*309.540 ##Jy
    
    return flux

def mom02SBgas(data):
    """
    data: incli(rad), redshift, bmaj, bmin, mom0
    from shi 2018, 1.36 include helium contribution
    return in Msun/pc^2
    
    """
    return 1.2*10**4*np.pi*np.cos(incli)*(1+redshift)**3*mom0/(bmaj*bmin)
    