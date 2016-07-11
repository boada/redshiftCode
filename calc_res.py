import numpy as np
import pyfits as pyf
from astLib import astSED
from scipy import optimize as opt

class Parameter:
    def __init__(self, value): self.value = value
    def set(self, value): self.value = value
    def __call__(self): return self.value

def fitter(function, parameters, y, x = None):
    def f(params):
        i=0
        for p in parameters:
            p.set(params[i])
            i+=1
        return y - function(x)
    if x is None: x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    opt.leastsq(f, p)

def gaussfit( x, y ):

    # GUESSING THE GAUSSIAN PARAMETERS
    mu0 = Parameter( (x[-1]+x[0])/2.0 )
    sig0 = Parameter( (x[-1]-x[0])/2.0 )
    h0 = Parameter( max(y)-min(y) )

    def gau(x): return h0() * np.exp(-((x-mu0())/sig0())**2)

    # PERFORMING FIT AND PRODUCING OUTPUTS
    fitter( gau, [mu0, sig0, h0], y, x )
    fitx = np.linspace( x[0], x[-1], 10*len(x) )
    fity = gau( fitx )

    return fitx, fity, [mu0.value, sig0.value, h0.value]


spectrum = 'bcs0000_imcmb3_oextr3_spectrum.fits'
#spectrum = 'bcs0324_imcmb1_oextr1_spectrum.fits'

data = pyf.getdata(spectrum)
oimg = pyf.open(spectrum)

owavelengthStep = oimg[1].header['CDELT1']
owavelengthStart = oimg[1].header['CRVAL1']

odata_dim = data.shape
owavelengthRange = np.array([owavelengthStart + i * owavelengthStep for i in
        range(odata_dim[0])])

linelist = np.loadtxt('linelist_may13')

total = 0
center = 0
res = 0
for line in linelist:
    mask = (line - 10 <owavelengthRange) & (owavelengthRange < line+10)
    x = owavelengthRange[mask]
    y = data[mask]

    _, _, params = gaussfit(x,y)

    print 'center', abs(line-params[0])
    print 'c res', abs(line-params[0]) * (3e5/line)
    print 'res', params[1] * (3e5/params[0])

    res += params[1] * (3e5/params[0])
    center += abs(line-params[0])
    total += np.sqrt((abs(line-params[0]) * (3e5/line))**2 + (params[1] *
    (3e5/params[0]))**2)


print '----'
print 'avg', total/len(linelist)
print 'center', center/len(linelist)
print 'res', res/len(linelist)
