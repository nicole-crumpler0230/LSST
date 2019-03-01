import sys, os
import math
import logging
import time
import yaml
import galsim as galsim
import galsim.wfirst as wfirst
import galsim.config.process as process
import galsim.des as des
from astropy.time import Time
import cProfile, pstats
import glob
import shutil
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
from astropy import units as u

def main(argv):
    pixel_scale = 0.2  # arcsec / pixel
    gal_image=galsim.Image(200,200,scale=pixel_scale)
    
    #Galaxy parameters
    bulge_n = 3.5          #
    bulge_re = 2.3         # arcsec
    disk_n = 1.5           #
    disk_r0 = 0.85         # arcsec
    bulge_frac = 0.3       #
    psf_sigma = 0.25       # arcsec

    # Define the galaxy profile
    bulge = galsim.Sersic(bulge_n, half_light_radius=bulge_re)
    disk = galsim.Sersic(disk_n, scale_radius=disk_r0)
    gal = bulge_frac * bulge + (1-bulge_frac) * disk
    
    # Define the PSF profile
    psf = galsim.Gaussian(flux=1., sigma=psf_sigma) # PSF flux should always = 1
    
    #Modifying the Galaxy with a given magnitude and redshift
    z=0.5 #redshift
    M=15 #magnitude
    filter_path= os.path.join(galsim.meta_data.share_dir,'bandpasses', 'LSST_r.dat')
    #gives path to location of bandpass data, different bandpass files can be found at
    #https://github.com/GalSim-developers/GalSim/tree/releases/2.1/share/bandpasses
    filters= galsim.Bandpass(filter_path, wave_type='nm') 
    filters=filters.withZeropoint('AB') #Set zeropoint of magnitude system
    sedpath_Gal   = os.path.join(galsim.meta_data.share_dir, 'SEDs', 'vega.txt')
    #gives path to location of SED data, different SED files can be found at
    #https://github.com/GalSim-developers/GalSim/tree/releases/2.1/share/SEDs
    sed   = galsim.SED(sedpath_Gal, wave_type='nm', flux_type='flambda')
    sed_ = sed.atRedshift(z)
    sed_ = sed_.withMagnitude(M, filters)
    gal=gal*sed_
    print (filters.effective_wavelength)
    
    #Make galaxy no longer have SED to create effective image at single wavelength
    flux = gal.calculateFlux(filters)
    gal= gal.evaluateAtWavelength(filters.effective_wavelength)
    gal= gal.withFlux(flux) # reapply correct flux
    
    #Convolve
    final = galsim.Convolve([gal, psf])
    
    #Create the image on a stamp
    b1=galsim.BoundsI(50,150,50,150)
    sub_gal_image_1=gal_image[b1]
    final.drawImage(sub_gal_image_1)

    if not os.path.isdir('Python_output'):
        os.mkdir('Python_output')
    file_name = os.path.join('Python_output', 'Galaxy_flux_redshift.fits')
    gal_image.write(file_name)
if __name__ == "__main__":
    main(sys.argv)
    import sys, os