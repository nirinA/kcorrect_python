'''Python wrapper for kcorrect
TODO:
. check ref count
. label goto fail
. option band_shift for filters
. option redshift
'''
import os
try:
    os.environ['KCORRECT_DIR']
except KeyError:
    raise SystemExit('Aborting! KCORRECT_DIR must be set')

import _kcorrect
import numpy

def load_templates(*args):
    vfile = "/mnt/big/store/kcorrect/data/templates/vmatrix.default.dat"
    lfile = "/mnt/big/store/kcorrect/data/templates/lambda.default.dat"
    _kcorrect.load_templates(vfile, lfile)

def load_filters(band_shift=0.):
    ffile = "/mnt/big/store/kcorrect/data/templates/sdss_filters.dat"
    _kcorrect.load_filters(ffile, band_shift)
   
def fit_coeffs_from_file(c, outfile = 'coeffs.dat'):
    return _kcorrect.fit_coeffs_from_file(c, outfile)

def fit_coeffs(c):
    c = numpy.array(c, dtype='float32')    
    if len(c) != 11:
        raise _kcorrect.error('incompatible number of arguments')
    return _kcorrect.fit_coeffs(c)

def reconstruct_maggies(c, all_redshift = -1.):
    c = numpy.array(c, dtype='float32')
    if len(c) != 6:
        raise _kcorrect.error('incompatible number of arguments')
    return _kcorrect.reconstruct_maggies(c, all_redshift)

def reconstruct_maggies_from_file(c, all_redshift = -1.):
    outfile = 'reconstructed_mag.dat'
    with open(outfile, 'w') as out:
        with open(c) as coeffs:
            for l in coeffs.readlines():
                l = list(map(float, l.strip('\n').split()))
                l = numpy.array(l, dtype='float32')
                r = reconstruct_maggies(l, all_redshift = -1.)
                out.write('%e %e %e %e %e %e\n'%tuple(r))

def fit_photoz(c):
    c = numpy.array(c, dtype='float32')
    return _kcorrect.fit_photoz(c)
