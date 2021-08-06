"""
pixel->sky;
sky->pixel;
read wcs from header
r&w ds9 region
"""

from photutils import EllipticalAperture
from astropy.io import fits
from astropy.wcs import WCS
from regions import PixCoord, EllipseSkyRegion, EllipsePixelRegion
from regions import read_ds9, write_ds9

def sum_aperture(x_c, y_c,a_out, b_out, PA,image):
    """
    Pa: degree
    """
    aperture=EllipticalAperture((x_c, y_c), a_out, b_out, theta=(PA+90)*np.pi/180) # mostly pixel
    phot_table = aperture_photometry(image, aperture_w1)
    phot= np.array(phot_table['aperture_sum'])
    
    return phot

def getwcs(filename):
    header=fits.getheader(filename)
    wcs=WCS(header)
    
    return wcs

def convert(type1,data,wcs):
    """
    p2s: pixel 2 sky;
    s3p: sky 2 pixel;
    data: xc,yc,bout(half),aout(half),pa (degree) for pixel input
    data: sky regions for sky input
    *only pixel data can be ploted*
    """
    if type1=='p2s':
        print('pixel to sky')
        ellipse_pix = EllipsePixelRegion(center=PixCoord(x=data[0], y=data[1]),
                                 height=data[2]*2, width=data[3]*2,angle=(data[4]+90) * u.deg)
        regions=ellipse_pix.to_sky(wcs=wcs)
        print('SkyCoord: ',regions)
    elif type1=='s2p':
        print('sky to pixel')
        pixel=data.to_pixel(wcs=wcs)
        # pixel.plot(ax,edgecolor)
        print('pixelCoord: ', pixel)
        
    return

def readds9(filename):
    regions = read_ds9(filename)
    return regions

def writeds9(regions,filename):
    write_ds9(regions, filename)
    return 

















        