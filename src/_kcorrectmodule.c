#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <kcorrect.h>

#include "numpy/npy_3kcompat.h"
#include "numpy/arrayobject.h"

static PyObject *_kcorrectError;

static IDL_LONG nz=1000;
static IDL_LONG maxiter=10000;
static IDL_LONG nprior=2;
static float tolerance=1.e-6;
static float zmin=0., zmax=1.e-0;

static float *lprior=NULL;
static float *zprior=NULL;
static float *lambda=NULL;
static float *vmatrix=NULL;
static float *rmatrix=NULL;
static float *zvals=NULL;
static IDL_LONG nk,nv,nl;
static float *filter_lambda=NULL;
static float *filter_pass=NULL;
static IDL_LONG *filter_n=NULL;
static IDL_LONG maxn,ndim,*sizes=NULL;
static float *redshift=NULL;
static float *maggies=NULL;
static float *maggies_ivar=NULL;
static float *coeffs=NULL;
static float *chi2=NULL;

float *pyvector_to_Carrayptrs(PyArrayObject *arrayin);

float *pyvector_to_Carrayptrs(PyArrayObject *arrayin)
{
    /* pointer to arrayin data as float */
    return (float *) arrayin->data; 
}

static PyObject *
kcorrect_load_templates(PyObject *self, PyObject *args)
{
    char * vfile;
    char * lfile;
    IDL_LONG k;
    if (!PyArg_ParseTuple(args, "ss", &vfile, &lfile))
        return NULL;

    /* read in templates */
    k = k_read_ascii_table(&vmatrix,&ndim,&sizes,vfile);
    if ((k!=1) && PyErr_Occurred()) {
	PyErr_SetString( _kcorrectError,"unable to read vmatrix template.\n");
	return NULL;
    }
    nl=sizes[1];
    nv=sizes[0];
    FREEVEC(sizes);
    k_read_ascii_table(&lambda,&ndim,&sizes,lfile);
    if ((sizes[0]!=nl+1) && PyErr_Occurred()) {
        PyErr_SetString( _kcorrectError,"incompatible vmatrix and lambda files.\n");
        return NULL;
    } /* end if */
    FREEVEC(sizes);

    Py_INCREF(Py_None);
    return Py_None;
}

PyDoc_STRVAR(kcorrect_load_templates_doc,
"loads the templates."
);

static PyObject *
kcorrect_load_filters(PyObject *self, PyObject *args)
{
    char * ffile;
    int i;
    float band_shift;
    if (!PyArg_ParseTuple(args, "sf", &ffile, &band_shift))
        return NULL;

    /* load in the filters */
    k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,ffile);

    /* create the rmatrix; this is a big matrix which tabulates the
 projection of each basis element b onto each filter k, as a
 function of redshift; you only have to project onto the filters
 here, since every other spectrum you will project onto the
 filters will be a linear combination of the basis elements b; you
 interpolate the rmatrix to get a particular redshift */
    rmatrix=(float *) malloc(nz*nv*nk*sizeof(float));
    zvals=(float *) malloc(nz*sizeof(float));
    for(i=0;i<nz;i++)
            zvals[i]=zmin+(zmax-zmin)*((float)i+0.5)/(float)nz;
    k_projection_table(rmatrix,nk,nv,vmatrix,lambda,nl,zvals,nz,filter_n,
                       filter_lambda,filter_pass,band_shift,maxn);
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyDoc_STRVAR(kcorrect_load_filters_doc,
"loads the filters."
);

static PyObject *
kcorrect_projection_table(PyObject *self, PyObject *args)
{
    PyArrayObject *py_lambda, *py_vmatrix, *py_rmatrix=NULL;
    npy_intp *shape;
    PyObject *newshape, *av, *al, *_iter = NULL;
    float *_lambda, *_vmatrix;
    float *_rmatrix=NULL;
    float *_zvals=NULL;
    IDL_LONG _nk,_nv,_nl;
    float *_filter_lambda=NULL;
    float *_filter_pass=NULL;
    IDL_LONG *_filter_n=NULL;
    IDL_LONG _maxn, _nz=1000;
    char * ffile;
    int i, ndim;
    float band_shift;
    npy_intp *dims=NULL;
    PyArray_Descr * dsc;
    if (!PyArg_ParseTuple(args, "O!O!fs|i",
			  &PyArray_Type, &py_vmatrix,
			  &PyArray_Type, &py_lambda,
			  &band_shift, &ffile, &_nz))
        return NULL;
    ndim = PyArray_NDIM(py_vmatrix);
    shape = PyArray_DIMS(py_vmatrix); /*PyArray_SHAPE(py_vmatrix); */
    _nv = shape[0];
    _nl = shape[1];
    newshape = Py_BuildValue("(i)", _nv*_nl);
    _lambda = pyvector_to_Carrayptrs(py_lambda);
    _vmatrix = pyvector_to_Carrayptrs((PyArrayObject *)PyArray_Reshape(py_vmatrix, newshape));
    k_load_filters(&_filter_n,&_filter_lambda,&_filter_pass,&_maxn,&_nk,ffile);
    _rmatrix=(float *) malloc(_nz*_nv*_nk*sizeof(float));
    dims=(npy_intp *) malloc(sizeof(npy_intp));    
    dims[0] = _nz*_nv*_nk;
    dsc = PyArray_DescrFromType(NPY_FLOAT);
    _zvals=(float *) malloc(_nz*sizeof(float));
    for(i=0;i<_nz;i++)
        _zvals[i]=zmin+(zmax-zmin)*((float)i+0.5)/(float)_nz;
    k_projection_table(_rmatrix,_nk,_nv,_vmatrix,_lambda,_nl,_zvals,_nz,_filter_n,
                       _filter_lambda,_filter_pass,band_shift,_maxn);

    py_rmatrix = (PyArrayObject *) PyArray_Zeros(1, dims, dsc, 0);
    if (NULL == py_rmatrix)  goto fail;
    _iter = PyArray_IterNew((PyObject*)py_rmatrix);
    if (_iter == NULL) goto fail;
    for(i=0;i<_nz*_nv*_nk;i++){
        if (PyArray_SETITEM(py_rmatrix, PyArray_ITER_DATA(_iter), PyFloat_FromDouble(_rmatrix[i]))) {
            PyErr_SetString(_kcorrectError, "unable to set _rmatrix");
            return NULL;
         }
         PyArray_ITER_NEXT(_iter);
    }

    FREEVEC(_rmatrix);
    FREEVEC(_zvals);
    FREEVEC(dims);
    Py_DECREF(_iter);
    Py_INCREF(py_rmatrix);
    return PyArray_Return(py_rmatrix);
fail:
    Py_XDECREF(_iter);
    Py_XDECREF(py_rmatrix);
    return NULL;


    /*return PyArray_Reshape(py_vmatrix, newshape);*/
    /*return PyArray_Return(py_rmatrix);*/
    /*return Py_BuildValue("ddi",_rmatrix[0], _vmatrix[1], dims[0]);*/
}

PyDoc_STRVAR(kcorrect_projection_table_doc,
"create the lookup table for given wavelength and flux."
);


static PyObject *
kcorrect_fit_coeffs_from_file(PyObject *self, PyObject *args)
{
    IDL_LONG i,j,k,niter,nchunk,ncurrchunk;
    float *mag;
    char * cfile, *ofile;
    PyArrayObject *cin = NULL;
    PyArray_Descr * dsc;
    FILE *fp, *ofp;
    int nd;

    if ((NULL == vmatrix) && (NULL == lambda)) {
            PyErr_SetString( _kcorrectError,"no templates loaded.\n");
            return NULL;;
    } /* end if */

    if (NULL == rmatrix) {
            PyErr_SetString( _kcorrectError,"no filters loaded.\n");
            return NULL;;
    } /* end if */

    if (!PyArg_ParseTuple(args, "ss", &cfile, &ofile))
        return NULL;
    fp=fopen(cfile,"r");
    ofp=fopen(ofile,"w");
    dsc = PyArray_DescrFromType(NPY_FLOAT32);
    cin = (PyArrayObject *)PyArray_FromFile(fp, dsc, -1, " ");
    fclose(fp);
    nd=cin->dimensions[0]/11;
    nchunk=nd;
    mag=pyvector_to_Carrayptrs(cin);
    redshift=(float *) malloc((nchunk+1)*sizeof(float));
    maggies=(float *) malloc((nchunk+1)*nk*sizeof(float));
    maggies_ivar=(float *) malloc((nchunk+1)*nk*sizeof(float));
    coeffs=(float *) malloc((nchunk+1)*nv*sizeof(float));
    chi2=(float *) malloc((nchunk+1)*sizeof(float));

    for (i=0;i<nd && i<nchunk;i++) {
        redshift[i] = mag[i*11];
        for(k=0;k<nk;k++) {
            maggies[i*nk+k] = mag[((i*11)+1)+k];
            maggies_ivar[i*nk+k] = mag[((i*11)+6)+k];
        }
    }

    ncurrchunk=nd;
    /* no direct constraints on the coeffs are included in this fit */
    k_fit_nonneg(coeffs,rmatrix,nk,nv,zvals,nz,maggies,
                 maggies_ivar,redshift,ncurrchunk,tolerance,
                 maxiter,&niter,chi2,0,0);
    for(i=0;i<ncurrchunk;i++) {
            fprintf(ofp,"%e ",redshift[i]);
            for(j=0;j<nv;j++)
                    fprintf(ofp,"%e ",coeffs[i*nv+j]);
            fprintf(ofp,"\n");
    } /* end for i */
    fclose(ofp);

    FREEVEC(redshift);
    FREEVEC(maggies);
    FREEVEC(maggies_ivar);
    FREEVEC(coeffs);
    FREEVEC(chi2);

    Py_INCREF(Py_None);
    return Py_None;
}

PyDoc_STRVAR(kcorrect_fit_coeffs_from_file_doc,
"writes the computed coeffs to a file\n\
with redshift and maggies read from a given file"
);

static PyObject *
kcorrect_fit_coeffs(PyObject *self, PyObject *args)
{
    IDL_LONG j,niter;
    PyArrayObject *pyin, *pyout;
    float *cin, *cout;
    npy_intp dims[] = {6};
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_FLOAT32);

    if ((NULL == vmatrix) && (NULL == lambda)) {
            PyErr_SetString( _kcorrectError,"no templates loaded.\n");
            return NULL;;
    } /* end if */

    if (NULL == rmatrix) {
            PyErr_SetString( _kcorrectError,"no filters loaded.\n");
            return NULL;;
    } /* end if */

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &pyin))
        return NULL;
    if (NULL == pyin)  return NULL;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,
                                                 dsc,
                                                 1,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cin=pyvector_to_Carrayptrs(pyin);
    cout=pyvector_to_Carrayptrs(pyout);

    redshift=(float *) malloc(sizeof(float));
    maggies=(float *) malloc(nk*sizeof(float));
    maggies_ivar=(float *) malloc(nk*sizeof(float));
    coeffs=(float *) malloc(nv*sizeof(float));
    chi2=(float *) malloc(sizeof(float));

    cout[0] = redshift[0] = cin[0];
    maggies[0] = cin[1];
    maggies[1] = cin[2];
    maggies[2] = cin[3];
    maggies[3] = cin[4];
    maggies[4] = cin[5];
    maggies_ivar[0] = cin[6];
    maggies_ivar[1] = cin[7];
    maggies_ivar[2] = cin[8];
    maggies_ivar[3] = cin[9];
    maggies_ivar[4] = cin[10];

    /* no direct constraints on the coeffs are included in this fit */
    k_fit_nonneg(coeffs,rmatrix,nk,nv,zvals,nz,
                 maggies,maggies_ivar,redshift,
                 1,tolerance,maxiter,&niter,chi2,0,0);

    for(j=0;j<nv;j++)   cout[1+j] = coeffs[j];

    FREEVEC(redshift);
    FREEVEC(maggies);
    FREEVEC(maggies_ivar);
    FREEVEC(coeffs);
    FREEVEC(chi2);
    
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

PyDoc_STRVAR(kcorrect_fit_coeffs_doc,
"computes the coeffs from given maggies and redshift"
);

static PyObject *
kcorrect_reconstruct_maggies(PyObject *self, PyObject *args)
{
    IDL_LONG j;
    PyArrayObject *pyin, *pyout;
    float *cin, *cout;
    float all_redshift;
    npy_intp dims[] = {6};
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_FLOAT32);

    if ((NULL == vmatrix) && (NULL == lambda)) {
            PyErr_SetString( _kcorrectError,"no templates loaded.\n");
            return NULL;;
    } /* end if */

    if (NULL == rmatrix) {
            PyErr_SetString( _kcorrectError,"no filters loaded.\n");
            return NULL;;
    } /* end if */

    if (!PyArg_ParseTuple(args, "O!f",
                          &PyArray_Type, &pyin,
                          &all_redshift))
        return NULL;
    if (NULL == pyin)  return NULL;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,             
                                                 dsc,
                                                 1,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cin=pyvector_to_Carrayptrs(pyin);
    cout=pyvector_to_Carrayptrs(pyout);

    redshift=(float *) malloc(sizeof(float));
    maggies=(float *) malloc(nk*sizeof(float));
    coeffs=(float *) malloc(nv*sizeof(float));
    if(all_redshift!=-1.) {
        redshift[0]=all_redshift;
    } else {
        redshift[0]=cin[0];
    }

    coeffs[0] = cin[1];
    coeffs[1] = cin[2];
    coeffs[2] = cin[3];
    coeffs[3] = cin[4];
    coeffs[4] = cin[5];
    k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,redshift,maggies,1);
    cout[0] = redshift[0];
    for(j=0;j<nv;j++)   cout[1+j] = maggies[j];

    FREEVEC(redshift);
    FREEVEC(maggies);
    FREEVEC(coeffs);
    
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

PyDoc_STRVAR(kcorrect_reconstruct_maggies_doc,
"reconstructs maggies from given coeffs and redshift"
);

static PyObject *
kcorrect_fit_photoz(PyObject *self, PyObject *args)
{
    IDL_LONG j,niter;
    PyArrayObject *pyin, *pyout;
    float *cin, *cout;
    npy_intp dims[] = {6};
    PyArray_Descr * dsc;
    dsc = PyArray_DescrFromType(NPY_FLOAT32);

    if ((NULL ==vmatrix) && (NULL == lambda)) {
            PyErr_SetString( _kcorrectError,"no templates loaded.\n");
            return NULL;;
    } /* end if */

    if (NULL == rmatrix) {
            PyErr_SetString( _kcorrectError,"no filters loaded.\n");
            return NULL;;
    } /* end if */

    if (!PyArg_ParseTuple(args, "O!",
                          &PyArray_Type, &pyin))
        return NULL;
    if (NULL == pyin)  return NULL;
    pyout=(PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type,             
                                                 dsc,
                                                 1,
                                                 dims,
                                                 NULL,
                                                 NULL,
                                                 0,
                                                 NULL);
    if (NULL == pyout)  return NULL;
    cin=pyvector_to_Carrayptrs(pyin);
    cout=pyvector_to_Carrayptrs(pyout);
    lprior=(float *) malloc(nprior*sizeof(float));
    zprior=(float *) malloc(nprior*sizeof(float));
    zprior[0]=0.;
    zprior[1]=1000.;
    redshift=(float *) malloc(sizeof(float));
    maggies=(float *) malloc(nk*sizeof(float));
    maggies_ivar=(float *) malloc(nk*sizeof(float));
    coeffs=(float *) malloc(nv*sizeof(float));
    chi2=(float *) malloc(sizeof(float));
    
    maggies[0] = cin[0];
    maggies[1] = cin[1];
    maggies[2] = cin[2];
    maggies[3] = cin[3];
    maggies[4] = cin[4];
    maggies_ivar[0] = cin[5];
    maggies_ivar[1] = cin[6];
    maggies_ivar[2] = cin[7];
    maggies_ivar[3] = cin[8];
    maggies_ivar[4] = cin[9];

    k_fit_photoz(redshift,coeffs,rmatrix,nk,nv,zprior,lprior,nprior,zvals, 
                 nz,maggies,maggies_ivar,1,tolerance,maxiter,&niter,
                 chi2,0);
    cout[0] = redshift[0];
    for(j=0;j<nv;j++)   cout[1+j] = coeffs[j];

    FREEVEC(redshift);
    FREEVEC(maggies);
    FREEVEC(maggies_ivar);
    FREEVEC(coeffs);
    FREEVEC(chi2);
    FREEVEC(lprior);
    FREEVEC(zprior);
    
    Py_INCREF(pyout); 
    return PyArray_Return(pyout);                     
}

PyDoc_STRVAR(kcorrect_fit_photoz_doc,
"returns guessed coeffs at redshift"
);

static PyObject *
kcorrect_k_template(PyObject *self, PyObject *args)
{
    float *v = NULL;
    char *vfile;
    PyArrayObject *pyout;
    PyArray_Descr *dsc;
    PyObject *pyout_iter = NULL;
    npy_intp *dims=NULL;
    IDL_LONG n, *s=NULL;
    int i,j=1;
    dsc = PyArray_DescrFromType(NPY_FLOAT);
    if (!PyArg_ParseTuple(args, "s", &vfile))
        return NULL;
    k_read_ascii_table(&v,&n,&s,vfile);
    dims = (npy_intp *) malloc((n)*sizeof(npy_intp));    
    if (NULL == dims) return NULL;
    for (i=0; i<n; i++) {
        dims[i] = s[i];
	j *= s[i];
    }
    pyout = (PyArrayObject *) PyArray_Zeros((npy_intp)n, dims, dsc, 0);
    if (NULL == pyout)  goto fail;
    pyout_iter = PyArray_IterNew((PyObject*)pyout);
    if (pyout_iter == NULL) {
        PyErr_SetString(_kcorrectError, "cannot create iterators");
        goto fail;
    }
    for (i=0;i<j;i++) {
      if (PyArray_SETITEM(pyout, PyArray_ITER_DATA(pyout_iter), PyFloat_FromDouble((double)v[i]))) {
	    PyErr_SetString(_kcorrectError, "unable to set vmatrix");
	    goto fail;
	}
        PyArray_ITER_NEXT(pyout_iter);
    }
    Py_DECREF(pyout_iter);
    Py_INCREF(pyout);
    return PyArray_Return(pyout);

fail:
    Py_XDECREF(pyout_iter);
    Py_XDECREF(pyout);
    return NULL;
}

PyDoc_STRVAR(kcorrect_k_template_doc,
"return vmatrix or lambda from a template file."
);

static PyMethodDef kcorrect_methods[] = {
    {"load_templates", kcorrect_load_templates, METH_VARARGS, kcorrect_load_templates_doc},
    {"load_filters", kcorrect_load_filters, METH_VARARGS, kcorrect_load_filters_doc},
    {"fit_coeffs_from_file", kcorrect_fit_coeffs_from_file, METH_VARARGS, kcorrect_fit_coeffs_from_file_doc},
    {"fit_coeffs", kcorrect_fit_coeffs, METH_VARARGS, kcorrect_fit_coeffs_doc},
    {"reconstruct_maggies", kcorrect_reconstruct_maggies, METH_VARARGS, kcorrect_reconstruct_maggies_doc},
    {"fit_photoz", kcorrect_fit_photoz, METH_VARARGS, kcorrect_fit_photoz_doc},
    {"k_template", kcorrect_k_template, METH_VARARGS, kcorrect_k_template_doc},
    {"projection_table", kcorrect_projection_table, METH_VARARGS, kcorrect_projection_table_doc},
    {NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module provides Python wrapper for kcorrect library"
);

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef kcorrectmodule = {
	PyModuleDef_HEAD_INIT,
	"kcorrect",
	module_doc,
	-1,
	kcorrect_methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC
PyInit__kcorrect(void)
{
	PyObject *m;
	m = PyModule_Create(&kcorrectmodule);  
	import_array();
        if (m == NULL)
            return NULL;
        _kcorrectError = PyErr_NewException("_kcorrect.error", NULL, NULL);
        Py_INCREF(_kcorrectError);
        PyModule_AddObject(m, "error", _kcorrectError);
        return m;
}
#else
PyMODINIT_FUNC
init_kcorrect(void)
{
	PyObject *m;
	m = Py_InitModule3("_kcorrect", kcorrect_methods, module_doc);
	import_array();
	if (m == NULL)
            goto finally;
        _kcorrectError = PyErr_NewException("_kcorrect.error", NULL, NULL);
        Py_INCREF(_kcorrectError);
        PyModule_AddObject(m, "error", _kcorrectError);
        finally:
        return;
}
#endif
