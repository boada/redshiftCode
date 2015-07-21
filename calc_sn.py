import pyfits as pyf
from astLib import astSED
from astLib import astWCS
from astLib import astStats
from astLib import astCoords
import numpy as np
import pylab as pl

def loadIFUSpectra(objectFileName):
    """ Loads in an object spectrum - this has to be in DEEP2 pipeline
    spec1d format (i.e. fits tables)
    Object spectrum is smoothed by boxcar of size smoothPix.

    Returns a dictionary containing object and sky astSED.SED objects
    {'object', 'sky'}

    """

    print "Loading IFU spectrum ..."

    oimg = pyf.open(objectFileName)

    objectFileName = objectFileName
    oimg = pyf.open(objectFileName)

    # Load the IFU data -- Row-stacked spectra
    odata = oimg[1].data
    oError = oimg[2].data
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
        #oflux = odata[i] - oskyflux
        oflux = odata[i]
        oErrorFlux = oError[i]

        # Mask out extreme values in spectrum
        # Just because edges dodgy in efosc
        med = np.median(oflux)
        oflux[np.greater(abs(oflux), 10.0*med)] = 0.0001

        objSED = astSED.SED(wavelength=owavelengthRange, flux=oflux)

        #  make it > 0 everywhere
        #objSED.flux = objSED.flux - objSED.flux.min()
        #objSED.flux = objSED.flux / objSED.flux.max()

        skySED = astSED.SED(wavelength=owavelengthRange, flux=oskyflux)
        errSED = astSED.SED(wavelength=owavelengthRange, flux=oErrorFlux)

        RSS.append({'object': objSED, 'sky': skySED, 'error' : errSED})

    return RSS

def maskLines(objSED, skySED):
    # Mask prominent sky emission lines
    if skySED is None:
        fixbad = "n"
    else:
        fixbad = "y"
    if fixbad == "y":
        normSkyFlux = skySED.flux / skySED.flux.max()
        threshold = 0.15
        badPix = np.where(normSkyFlux > threshold)[0]
        lines = []
        for i in range(len(badPix)):
            startPixInLine = False
            for line in lines:
                if line[0] <= badPix[i] <= line[1]:
                    startPixInLine = True
            if not startPixInLine:
                pixCount = 1
                if pixCount + i < len(badPix) - 1:
                    nextPix = badPix[i + pixCount]
                    prevPix = badPix[i]
                    maxReached = False
                    while nextPix < prevPix + 2:
                        if pixCount + i < len(badPix) - 1:
                            prevPix = badPix[i + pixCount]
                            pixCount += 1
                            nextPix = badPix[i + pixCount]
                        else:
                            maxReached = True
                            break
                    if not maxReached:
                        lastPix = prevPix
                    else:
                        lastPix = max(badPix)
                else:
                    lastPix = max(badPix)
                #append the lines and add a little padding
                lines.append([badPix[i]-5, lastPix+5])

        for line in lines:
            # Do a simple linear fit to the end points
            y = objSED.flux[line[1]] - objSED.flux[line[0]]
            x = skySED.wavelength[line[1]] - skySED.wavelength[line[0]]
            slope = y/x
            intercept = objSED.flux[line[0]] - slope *\
                skySED.wavelength[line[0]]
#                print skySED.wavelength[line[0]], skySED.wavelength[line[1]]

            for i in range(line[0], line[1]):
                objSED.flux[i] = slope * skySED.wavelength[i] + intercept

        return objSED

def median_continuum(flux, error, numsig=1.5, plot=False):
    """ Estimates the continuum using a median and sigma clipping.

    Given the fluxes and one sigma errors for the section, calculates
    the flux median. Then rejects all flux values less than numsig*sig
    lower than the median, where sig is the median one sigma error
    value of the flux values. Repeat this process, only retaining the
    not-rejected flux values each time, until no flux values are
    rejected. Once this condition is reached, take the current median
    value as the continuum.

    Returns the continuum value.
    """

    if plot:
        pl.plot(flux)
    while True:
        medfl = np.median(flux)
        meder = np.median(error)
        if plot:
            l = pl.axhline(medfl)
        cond = (flux > (medfl - meder*numsig))
        badflux = flux[~cond]
        if len(badflux) == 0:
            return medfl
        flux = flux[cond]
        error = error[cond]


objectFileName = 'bcs0267_crcl_oextr3.fits'

RSS = loadIFUSpectra(objectFileName)


obj = RSS[246-79]['object']
sky = RSS[246-79]['sky']
err = RSS[246-79]['error']

#obj = maskLines(obj,sky)
#sky = maskLines(sky, sky)


sn = sum(obj.flux)/np.sqrt(sum(sky.flux))
print sn


