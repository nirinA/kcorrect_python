'''some utilities to manipulate magnitudes...

solar_magnitude is adapted from Python code originally
written by Benjamin Weaver.
'''
import os
import numpy
import kcorrect

## fit_nonneg expects arguments in maggy units
## may also need nanomaggy converter...
def mag2maggies(mag):
    return numpy.power(10., -0.4*mag)

def invariance(maggies, mag_err):
    return  numpy.power(0.4*numpy.log(10.)*maggies*mag_err, -2)

def read_basel(solarname, silent=False):
    """
    Read a spectrum from a Basel spectrum file.
    """
    filename = os.path.join(
        os.getenv('KCORRECT_DIR'), 'data/basel', solarname)
    if not os.path.isfile(filename):
        raise _kcorrect.error('%s not found'%filename)
    with open(filename,'r') as f:
        returndata = {'flux':[], 'model':[], 'teff':[], 'logg':[], 'mh':[],
            'vturb':[], 'xh':[]}
        rawdata = f.read()
        alldata = rawdata.replace("\n", '').split()
        returndata['wavelength'] = numpy.array(alldata[0:1221], dtype=numpy.float32)
        del alldata[0:1221]
        nunits = int(len(alldata)/(1227))
        if not silent:
            print ("%d block(s) of spectra" % nunits)
        for u in range(nunits):
            returndata['model'].append(int(alldata[0]))
            returndata['teff'].append(int(alldata[1]))
            returndata['logg'].append(float(alldata[2]))
            returndata['mh'].append(float(alldata[3]))
            returndata['vturb'].append(float(alldata[4]))
            returndata['xh'].append(float(alldata[5]))
            returndata['flux'].append(numpy.array(alldata[6:1227], dtype=numpy.float32))
            del alldata[0:1227]
    return returndata

def wavelength_to_edges(centers):
    """
    Convert set of pixel centers to equivalent edges.
    """
    n = len(centers)
    edges = numpy.zeros((n+1,),centers.dtype)
    edges[1:n] = 0.5*(centers[0:(n-1)] + centers[1:n])
    edges[0] = centers[0] - (edges[1] - centers[0])
    edges[n] = centers[n-1] + (centers[n-1] - edges[n-1])
    return edges
    
def solar_magnitudes(solarname='lcbsun.ori',
                     band_shift= 0.0,
                     filterlist='sdss_filters.dat',
                     silent=False):
    """
    Calculate the Solar absolute magnitudes in various bandpasses.
    """
    #
    # Read in the Sun & put it at 10 pc.
    #
    basel = read_basel(solarname)
    nspectra = len(basel['model'])
    angstroms = basel['wavelength']*10.0
    c = 2.99792458e+18 # angstrom sec^-1
    pctocm = 3.085677581e+18 # cm pc^-1
    Rsun = 6.96e+10 # cm
    flux = numpy.array([ 4.0*numpy.pi*f*c*((0.1*Rsun/pctocm)**2)/(angstroms**2) for f in numpy.array(basel['flux']) ], dtype='float32')
    edges = wavelength_to_edges(angstroms)
    maggies = kcorrect.project_filters(edges, flux, band_shift=band_shift, filterlist=filterlist)
    mag = -2.5*numpy.log10(maggies.reshape((maggies.size,)))
    return mag
