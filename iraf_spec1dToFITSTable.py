#!/usr/bin/env python

"""Convert IRAF spec1d format to FITS tables used by my redshift code.

"""

import os
import sys
import numpy
from numpy import fft
import pylab
import pyfits
import glob
from astLib import *
from scipy import interpolate
from scipy import ndimage
from scipy import stats
from scipy import signal
from pyraf import iraf
from iraf import pysalt
from iraf import saltspec
from iraf import twodspec
from iraf import longslit
from iraf import onedspec
from iraf import apextract
import IPython
pylab.matplotlib.interactive(True)

#------------------------------------------------------------------------------------------------------------------
def onedspecToBintables(inFileName, outFileName, skyRow = 2):
    """Converts a file in iraf onedspec format into bintables format (for use with redshift code)
    
    skyRow gives the row in which the sky appears (row 1 is for variance in current apall set up)
    
    """

    img=pyfits.open(inFileName)
    data=img[0].data
    header=img[0].header
    refWavelength=header['CRVAL1']
    refPix=header['CRPIX1']
    angsPerPix=header['CD1_1']
    if data.ndim == 2:
        angstroms=(numpy.arange(0, data[0].shape[1], dtype=float)-refPix)*angsPerPix+refWavelength
    elif data.ndim == 1:
        angstroms=(numpy.arange(0, data.shape[0], dtype=float)-refPix)*angsPerPix+refWavelength
    else:
        raise Exception, "3d spectrum? We didn't decide how to treat those yet."
    
    # Put blue end at right hand side
    if angstroms[1] < angstroms[0]:
        angstroms=angstroms[::-1]
        data=data[::-1]
                
    # Save as .fits table - here we're assuming we have sky spectra in data[1], object in data[0]
    if data.ndim == 2:
        specColumn=pyfits.Column(name='SPEC', format='D', array=data[0][0])
        skyColumn=pyfits.Column(name='SKYSPEC', format='D', array=data[skyRow][0])
    elif data.ndim == 1:
        specColumn=pyfits.Column(name='SPEC', format='D', array=data)
        skyColumn=pyfits.Column(name='SKYSPEC', format='D', array=numpy.zeros(data.shape))
    lambdaColumn=pyfits.Column(name='LAMBDA', format='D', array=angstroms)
    tabHDU=pyfits.new_table([specColumn, lambdaColumn, skyColumn])
    tabHDU.name='1D_SPECTRUM'
    HDUList=pyfits.HDUList([pyfits.PrimaryHDU(), tabHDU])
    HDUList.writeto(outFileName, clobber=True)

#------------------------------------------------------------------------------------------------------------------    
# Main
if len(sys.argv) < 3:
    print "Run: % iraf_spec1dToFITSTable.py <input iraf spec1d .fits> <output FITS table .fits>"
else:
    
    inFileName=sys.argv[1]
    outFileName=sys.argv[2]
    
    onedspecToBintables(inFileName, outFileName)
    
    

    