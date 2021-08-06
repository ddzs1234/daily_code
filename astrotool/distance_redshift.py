from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.cosmology import  z_at_value
from astropy.coordinates import Distance

def distance2red(distance,zmax,ho,Om0):
    '''
	distance in Mpc;
	zmax
	h0: like 70,73,100
	Om0: matter 
	flat
    '''

    # h=0.73, Omega_0=0.27, omega_lambd=0.73 from shi 2021
    cosmo = FlatLambdaCDM(H0=h0, Om0=Om0, Tcmb0=2.725)
    print(cosmo)
    a = z_at_value(cosmo.luminosity_distance, distance * u.Mpc,zmax=zmax)

    return a 

def red2distance(z):
    """
    https://docs.astropy.org/en/stable/api/astropy.coordinates.Distance.html
    """
    return Distance(z)