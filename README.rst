what is kcorrect_python?
========================

This package provides Python interface to kcorrect C library by 
M. Blanton et al. and described here:

 `http://adsabs.harvard.edu/abs/2007AJ....133..734B<http://adsabs.harvard.edu/abs/2007AJ....133..734B>`_

requirements
============

kcorrect
--------

see:
    
  http://howdy.physics.nyu.edu/index.php/Kcorrect

for obtaining kcorrect and how to install it.

This wrapper uses kcorrect version 4.2.
  
Python and dependencies
-----------------------

This version requires Python 2.7+ and NumPy 1.7+

installation
============

the usual::

    python setup.py build

and (may need root privileges) ::

    python setup.py install

should build the package and install *_kcorret.so* and *kcorrect/*
into the standard *site-packages* directory.

usage
=====

Note: the environmental variables *KCORRECT_DIR* and
*LD_LIBRARY_PATH* should be set and point to the location
where the kcorrect package is installed, eg::

    export KCORRECT_DIR=/usr/local/kcorrect
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$KCORRECT_DIR/lib

test
----

you can test the package by running::

    python test.py

this will procude two files *coeffs.dat* and
*reconstructed_maggies.dat".

Available functions
-------------------

The following functions are currently available for this version:

    o :func:`load_templates` 
    o :func:`load_filters` 
    o :func:`fit_coeffs_from_file` 
    o :func:`fit_coeffs` 
    o :func:`reconstruct_maggies` 
    o :func:`reconstruct_maggies_from_files` 
    o :func:`fit_photoz` 
    o :func:`fit_photoz_from_file`
    o :func:`solar_magnitudes`
    o :func:`project_filters`
    o :func:`fit_nonneg`


LICENSE
=======

Public Domain

