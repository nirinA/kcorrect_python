'''Python wrapper for kcorrect

example usage

03:04:21$ python3 -c "import kcorrect; kcorrect.load_templates();kcorrect.load_filters();import numpy; a=numpy.array([0.03077382, 1.144068e-08, 5.262234e-08, 8.210213e-08, 8.744532e-08, 1.017738e-07, 6.216309e+16, 3.454767e+17, 1.827409e+17, 1.080889e+16, 3163927000000000.0],dtype='float32');print(a);b=kcorrect.fit_coeffs(a);print(b)"
[  3.07738204e-02   1.14406804e-08   5.26223403e-08   8.21021331e-08
   8.74453221e-08   1.01773800e-07   6.21630889e+16   3.45476690e+17
   1.82740894e+17   1.08088897e+16   3.16392712e+15]
[  3.07738204e-02   2.02254747e-14   1.49129165e-35   2.15513887e-06
   6.94462278e-06   1.78061924e-13]

>>> 
'''
import os
try:
    os.environ['KCORRECT_DIR']
except KeyError:
    raise SystemExit('Aborting! KCORRECT_DIR must be set')

import _kcorrect
import numpy as np

def load_templates(*args):
    vfile = "/mnt/big/store/kcorrect/data/templates/vmatrix.default.dat"
    lfile = "/mnt/big/store/kcorrect/data/templates/lambda.default.dat"
    _kcorrect.load_templates(vfile, lfile)

def load_filters(*args):
    ffile = "/mnt/big/store/kcorrect/data/templates/sdss_filters.dat"
    _kcorrect.load_filters(ffile)
   
def fit_coeffs_from_file(c, *args):
    outfile = 'coeffs.dat'
    return _kcorrect.fit_coeffs_from_file(c, outfile)

def fit_coeffs(c):
    return _kcorrect.fit_coeffs(np.array(c, dtype='float32'))
