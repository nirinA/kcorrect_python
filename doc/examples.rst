Examples of usage
-----------------

Note: The examples below use the data shipped with kcorrect.v4_2.

fit_coeffs and fit_photoz
^^^^^^^^^^^^^^^^^^^^^^^^^

You can use these functions as follow::

    >>> import kcorrect, numpy
    >>> kcorrect.load_templates()
    >>> kcorrect.load_filters()
    >>> a = numpy.array([0.03077382, 1.144068e-08, 5.262234e-08, 8.210213e-08, 8.744532e-08, 1.017738e-07, 6.216309e+16, 3.454767e+17, 1.827409e+17, 1.080889e+16, 3163927000000000.0], dtype='float32')
    >>> c = kcorrect.fit_coeffs(a)
    >>> c
    array([  3.07738204e-02,   2.02254747e-14,   1.49129165e-35,
         2.15513887e-06,   6.94462278e-06,   1.78061924e-13], dtype=float32)
    >>> m = kcorrect.reconstruct_maggies(c)
    >>> m
    array([  3.07738204e-02,   1.44426586e-08,   5.28384980e-08,
         8.09117182e-08,   9.51680121e-08,   1.10408600e-07], dtype=float32)

The example above successively loads the module,
loads the default templates, *vmatrix.default.dat*
and *lambda.default.dat*, loads the default filter,
*sdss_filters.dat*, then computes the coeffs and
reconstructs maggies.

The argument, ``a``, of :func:`fit_coeffs` consists of::

    >>> a[0] # redshift
    >>> a[1:6] # maggies
    >>> a[6:12] # maggies_ivar

To compute the reconstructed maggies at rest-frame with bandpasses
shifted by 0.1, you need first reload the filters with the given
band_shift, then compute the coeffs and the maggies::
    
    >>> kcorrect.load_filters(band_shift=0.1)
    >>> m0 = kcorrect.reconstruct_maggies(c, redshift=0.)

If the redshifs, maggies and maggies_invvar are stored
in a file like *sample.dat* found in the *test* directory
of kcorrect package, you can use :func:`fit_coeffs_from_file`
and :func:`reconstruct_maggies_from_files`  to perform the
computation::

    >>> kcorrect.fit_coeffs_from_file('some_file.dat', outfile='output_coeffs.dat')
    >>> kcorrect.reconstruct_maggies_from_files('output_coeffs.dat', outfile='computed_maggies.dat')

these produce 2 files *output_coeffs.dat* and *computed_maggies.dat*

To use different templates, you load them as follow::
    
    >>> kcorrect.load_templates(v='vmatrix.goods.dat',l='lambda.goods.dat')

If templates and filters are not loaded before calling the other
functions, error is raised::
    
    >>> kcorrect.fit_coeffs(range(11))
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "kcorrect.py", line 37, in fit_coeffs
        return _kcorrect.fit_coeffs(c)
    _kcorrect.error: no filters loaded.

:func:`fit_photoz` and  :func:`fit_photoz_from_file` can be used
as follow, after loading the appropriate templates and filter::

    >>> p = kcorrect.fit_photoz(a[1:])
    >>> p
    array([  1.41886109e-02,   5.18920551e-09,   6.65258128e-36,
         2.18073205e-06,   5.97664302e-06,   4.88666385e-14], dtype=float32)

if the data are from a file *photoz.dat*::

    >>> fit_photoz_from_file('photoz.dat', outfile='photoz.out')

which produces the result to the output file *photoz.out*

computing the k-correction
^^^^^^^^^^^^^^^^^^^^^^^^^^

first, load templates and filters, with defaults arguments here::
    
    >>> kcorrect.load_templates()
    >>> kcorrect.load_filters()

then, with the following inputs, compute the coeffs::
    
    >>> redshift = 0.03077382
    >>> maggies = numpy.array([1.144068e-08, 5.262234e-08, 8.210213e-08, 8.744532e-08, 1.017738e-07], dtype='float32')
    >>> maggies_ivar = numpy.array([6.216309e+16, 3.454767e+17, 1.827409e+17, 1.080889e+16, 3163927000000000.0], dtype='float32')
    >>> coeffs = kcorrect.fit_nonneg(redshift, maggies, maggies_ivar)
    >>> coeffs
    array([  3.07738204e-02,   2.02254747e-14,   1.49129165e-35,
             2.15513887e-06,   6.94462278e-06,   1.78061924e-13], dtype=float32)

the reconstructed maggies::

    >>> rm = kcorrect.reconstruct_maggies(coeffs)
    >>> rm
    array([  3.07738204e-02,   1.44426586e-08,   5.28384980e-08,
             8.09117182e-08,   9.51680121e-08,   1.10408600e-07], dtype=float32)

reload the filters with the appropriate *band_shift*, if needed::

    >>> kcorrect.load_filters(band_shift=0.1)

finally, compute the reconstructed maggies at redshit 0.::

    >>> rm0 = kcorrect.reconstruct_maggies(coeffs, redshift=0.)

the kcorrection is then::

    >>> kc = -2.5*numpy.log10(rm[1:]/rm0[1:])
    >>> kc
    array([ -3.17973822e-01,  -2.27918401e-01,  -1.61477402e-01,
             1.88272024e-04,  -9.79175493e-02], dtype=float32)
