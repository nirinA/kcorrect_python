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
    return _kcorrect.fit_photoz(c)

def fit_photoz_from_file(c, outfile='photoz_coeffs.dat'):
    with open(outfile, 'w') as out:
        with open(c) as maggies:
            for l in maggies.readlines():
                l = list(map(float, l.strip('\n').split()))
                l = numpy.array(l, dtype='float32')
                r = fit_photoz(l)
                out.write('%e %e %e %e %e %e\n'%tuple(r))

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
