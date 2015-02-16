import pyfits as pyf
from astLib import astSED
from astLib import astWCS
from astLib import astStats
from astLib import astCoords
import numpy as np

def loadIFUSpectra(objectFileName):
    """ Loads in an object spectrum - this has to be in DEEP2 pipeline
    spec1d format (i.e. fits tables)
    Object spectrum is smoothed by boxcar of size smoothPix.

    Returns a dictionary containing object and sky astSED.SED objects
    {'object', 'sky'}

    """

    print "Loading IFU spectrum ..."

    oimg = pyf.open(objectFileName)

    # Load the IFU data -- Row-stacked spectra
    odata = oimg[1].data
    odata_dim = odata.shape
    wcs = astWCS.WCS(objectFileName, extensionName=1)
    owavelengthStartEnd = wcs.getImageMinMaxWCSCoords()[0:2]
    fiberNumber = wcs.getImageMinMaxWCSCoords()[2:4]
    owavelengthStep = oimg[1].header['CDELT1']

    owavelengthRange = [owavelengthStartEnd[0] + i * owavelengthStep
                        for i in range(odata_dim[1])]

    # Check to make sure we got it right
    if not owavelengthRange[-1] == owavelengthStartEnd[-1]:
        print 'The ending wavelenghts do not match... Exiting'
        sys.exit(1)
    else:
        sums = [sum(odata[i,:]) for i in range(odata.shape[0])]
        #find the median value of all the fibers
        med = astStats.clippedMedianStdev(sums)
        med = med['clippedMedian']

        skyfibers = [i for i in range(odata.shape[0])\
                if sum(odata[i,:]) <= med]
        skydata = odata.take(skyfibers, axis=0)

        oskyflux = [np.average(skydata[:,i])\
                for i in range(skydata.shape[1])]

    RSS = []
    for i in range(int(fiberNumber[1])):
        oflux = odata[i] - oskyflux

        # Mask out extreme values in spectrum
        # Just because edges dodgy in efosc
        med = np.median(oflux)
        oflux[np.greater(abs(oflux), 10.0*med)] = 0.0001

        objSED = astSED.SED(wavelength=owavelengthRange, flux=oflux)

        #  make it > 0 everywhere
        #objSED.flux = objSED.flux - objSED.flux.min()
        #objSED.flux = objSED.flux / objSED.flux.max()

        skySED = astSED.SED(wavelength=owavelengthRange, flux=oskyflux)

        RSS.append({'object': objSED, 'sky': skySED})
    return RSS

#points = [(203.86237,41.013417),
#        (203.86496,41.013389),
#        (203.86104,41.011667),
#        (203.86375,41.011694),
#        (203.86633,41.011694),
#        (203.86225,41.009944),
#        (203.865,41.009944)]
#star = (203.861012858877,41.0105477333859)
#fibers = 129, 130, 143, 144, 145, 158, 159

points = [(203.81879,40.989833),
        (203.82137,40.989917),
        (203.82004,40.988222)]
star = 203.81698, 40.98942
fibers = 97, 98, 112

#points = [(203.80817,41.023667), (203.80679,41.021972), (203.8095,41.021944)]
#star = 203.806276951596, 41.0219357656736
#fibers = 35, 49, 50

objectFileName = 'bcs0045_crcl_oextr1.fits'

RSS = loadIFUSpectra(objectFileName)

x = [p[0] for p in points]
y = [p[1] for p in points]

w = []
for f in fibers:
    d = [flux for flux, wl in zip(RSS[246-f]['object'].flux,
        RSS[246-f]['object'].wavelength) if 5000<wl<5010]
    w.append(sum(d))

c = (np.average(x, weights=w), np.average(y, weights=w))
print w
print c

Ra_shift = astCoords.calcAngSepDeg(c[0],c[1],star[0],c[1])*3600
Dec_shift = astCoords.calcAngSepDeg(c[0],c[1],c[0],star[1])*3600

# Need to shift to the South West, so make the shifts negative
new_coords = astCoords.shiftRADec(c[0], c[1], -1*Ra_shift, -1*Dec_shift)
print Ra_shift, Dec_shift
print new_coords


