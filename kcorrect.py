'''Python wrapper for kcorrect

example usage

>>> import kcorrect
>>> kcorrect.fit_coeffs('/mnt/big/store/kcorrect_python/maggies.dat')
1.000000e-01 4.203895e-45 1.653600e+01 6.726233e-44 8.407791e-45 1.821688e-44 
1.000000e+00 3.082857e-44 3.704808e-01 8.371582e-38 4.344025e-44 7.103386e+01 
1.100000e+01 4.764415e-44 2.280125e+00 2.349515e-26 1.098618e-42 1.082857e+02 
>>> 
'''

import _kcorrect
import numpy as np

def read_maggies(mag):
    redshifts = []
    maggies = []
    maggies_ivar = []
    with open(mag) as dat:
        for l in dat.readlines():
            l = l.rstrip('\n').split()
            r,m,m_i = l[0], l[1:6], l[6:]
            redshifts.append(float(r))
            maggies.append(list(map(float,m)))
            maggies_ivar.append(list(map(float,m_i)))
    return {'redshifts':redshifts,
            'maggies':maggies,
            'maggies_ivar':maggies_ivar}

def load_templates(*args):
    vfile = "/mnt/big/store/kcorrect/data/templates/vmatrix.default.dat"
    lfile = "/mnt/big/store/kcorrect/data/templates/lambda.default.dat"
    _kcorrect.load_templates(vfile, lfile)

def load_filters(*args):
    ffile = "/mnt/big/store/kcorrect/data/templates/sdss_filters.dat"
    _kcorrect.load_filters(ffile)
   
def fit_coeffs(c, *args):
##    rmm = read_maffies(cfile)
##    redshifts = rmm['redshifts']
##    maggies = rmm['maggies']
##    maggies_ivar = rmm['maggies_ivar']
    return _kcorrect.fit_coeffs(c)
