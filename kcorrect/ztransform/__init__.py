import _ztransform
__doc__ = _ztransform.__doc__

__all__ = ['distmod']

def distmod(z,  omega0=0.3, omegal0=0.7):
    '''Calculate distance modulus."
INPUTS:
    z           redshifts
OPTIONAL INPUTS:
    omega0      omega_matter to use (default: 0.3)
    omegal0     omega_lambda to use (default: 0.7)
OUTPUTS:
    dm          distance modulus'''
    return _ztransform.z2dm(z, omega0, omegal0)

    
