'''Python wrapper for kcorrect
'''

import _kcorrect
import numpy as np


def fit_coeffs(c, *args):
    return _kcorrect.fit_coeffs(c)
