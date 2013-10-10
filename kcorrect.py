'''Python wrapper for kcorrect
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
    band_shift = 0
    _kcorrect.load_templates(vfile, lfile)

def load_filters(*args):
    ffile = "/mnt/big/store/kcorrect/data/templates/sdss_filters.dat"
    _kcorrect.load_filters(ffile)
   
def fit_coeffs_from_file(c, *args):
    outfile = 'coeffs.dat'
    return _kcorrect.fit_coeffs_from_file(c, outfile)

def fit_coeffs(c):
    return _kcorrect.fit_coeffs(numpy.array(c, dtype='float32'))

def reconstruct_maggies(c):
    band_shift = 0.
    redshift = 0.
    c = numpy.array(c, dtype='float32')
    return _kcorrect.reconstruct_maggies(c) ##, band_shift, redshift)

def reconstruct_maggies_from_file(c):
    outfile = 'reconstructed_mag.dat'
    with open(outfile, 'w') as out:
        with open(c) as coeffs:
            for l in coeffs.readlines():
                l = list(map(float, l.strip('\n').split()))
                l = numpy.array(l, dtype='float32')
                r = reconstruct_maggies(l)
                out.write('%e %e %e %e %e %e\n'%tuple(r))

def fit_photoz(c):
    c = numpy.array(c, dtype='float32')
    return _kcorrect.fit_photoz(c)
