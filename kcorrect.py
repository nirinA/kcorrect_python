'''Python wrapper for kcorrect
'''
import os
try:
    os.environ['KCORRECT_DIR']
except KeyError:
    raise SystemExit('Aborting! KCORRECT_DIR must be set')

import _kcorrect
import numpy

def load_templates(v="vmatrix.default.dat",
                   l="lambda.default.dat",
                   templates_dir='data/templates'):
    vfile = os.path.join(os.environ['KCORRECT_DIR'], templates_dir, v)
    if not os.path.isfile(vfile):
        raise _kcorrect.error('file: %s not found'%vfile)
    lfile = os.path.join(os.environ['KCORRECT_DIR'], templates_dir, l)
    if not os.path.isfile(lfile):
        raise _kcorrect.error('file: %s not found'%lfile)
    _kcorrect.load_templates(vfile, lfile)

def load_filters(f="sdss_filters.dat", band_shift=0.,
                 filters_dir='data/templates'):
    ffile = os.path.join(os.environ['KCORRECT_DIR'], filters_dir, f)
    if not os.path.isfile(ffile):
        raise _kcorrect.error('file: %s not found'%ffile)
    _kcorrect.load_filters(ffile, band_shift)
   
def fit_coeffs_from_file(c, outfile = 'coeffs.dat'):
    return _kcorrect.fit_coeffs_from_file(c, outfile)

def fit_coeffs(c):
    c = numpy.array(c, dtype='float32')    
    if len(c) == 11:
        return _kcorrect.fit_coeffs(c)
    else:
        raise _kcorrect.error('for this deault sdss filters, number of argument should be 11')

def reconstruct_maggies(c, redshift=-1.):
    c = numpy.array(c, dtype='float32')
    return _kcorrect.reconstruct_maggies(c, redshift)

def reconstruct_maggies_from_file(c, redshift=-1., outfile='reconstructed_mag.dat'):
    with open(outfile, 'w') as out:
        with open(c) as coeffs:
            for l in coeffs.readlines():
                l = list(map(float, l.strip('\n').split()))
                l = numpy.array(l, dtype='float32')
                r = reconstruct_maggies(l, redshift)
                out.write(' '.join([str(i) for i in r]))
                out.write('\n')

def fit_photoz(c):
    c = numpy.array(c, dtype='float32')
    if len(c) != 10:
        raise _kcorrect.error('incompatible number of arguments! should be 10')
    return _kcorrect.fit_photoz(c)

def fit_photoz_from_file(c, outfile='photoz_coeffs.dat'):
    with open(outfile, 'w') as out:
        with open(c) as maggies:
            for l in maggies.readlines():
                l = list(map(float, l.strip('\n').split()))
                l = numpy.array(l, dtype='float32')
                r = fit_photoz(l)
                out.write('%e %e %e %e %e %e\n'%tuple(r))

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

def project_filters(wavelength, flux, band_shift=0.,
                     filterlist='sdss_filters.dat',
                     filterlist_dir='data/templates',
                     zmin=0., zmax=1.,nz=100):
    filename = os.path.join(
        os.getenv('KCORRECT_DIR'), filterlist_dir, filterlist)
    if not os.path.isfile(filename):
        raise _kcorrect.error('filter list %s not found.'%filename)
    return _kcorrect.projection_table(wavelength, flux, \
                                      band_shift, \
                                      filename, \
                                      zmin, zmax, nz)
    
def solar_magnitudes(solarname='lcbsun.ori',
                     band_shift= 0.0,
                     filterlist='sdss_filters.dat',silent=False):
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
    maggies = project_filters(edges, flux, band_shift=band_shift, filterlist=filterlist)
    mag = -2.5*numpy.log10(maggies.reshape((maggies.size,)))
    return mag

def fit_nonneg(redshift, maggies, maggies_ivar):
    '''Fit nonnegative sum of given templates to given set of maggies'''
    redshift = numpy.array(redshift, dtype=numpy.float32)
    maggies = numpy.array(maggies, dtype=numpy.float32)
    maggies_ivar = numpy.array(maggies_ivar, dtype=numpy.float32)
    return _kcorrect.fit_nonneg(redshift, maggies, maggies_ivar)

def template(t, templates_dir='data/templates'):
    tfile = os.path.join(os.environ['KCORRECT_DIR'], templates_dir, t)
    if not os.path.isfile(tfile):
        raise _kcorrect.error('file: %s not found'%tfile)
    return _kcorrect.template(tfile)

    
