#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <kcorrect.h>

static PyObject *kcorrectError;

static IDL_LONG nz=1000;
static IDL_LONG maxiter=10000;
static float tolerance=1.e-6;
static float zmin=0., zmax=1.e-0;
static float band_shift=0.;

static float *lambda=NULL;
static float *vmatrix=NULL;
static float *rmatrix=NULL;
static float *zvals=NULL;
static IDL_LONG nk,nv,nl;
static float *filter_lambda=NULL;
static float *filter_pass=NULL;
static IDL_LONG *filter_n=NULL;
static IDL_LONG maxn;
static float *redshift=NULL;
static float *maggies=NULL;
static float *maggies_ivar=NULL;
static float *coeffs=NULL;
static float *chi2=NULL;

static PyObject *
kcorrect_fit_coeffs(PyObject *self, PyObject *args)
{
    IDL_LONG i,j,k,c,ndim,niter,nchunk,ncurrchunk,*sizes=NULL;
    const char * vfile = "vmatrix.default.dat";
    const char * lfile = "lambda.default.dat";
    const char * ffile = "sdss_filters.dat";
    const char * cfile = "";
    char * path = getenv("KCORRECT_DIR");
    const char * vmatrixfile;
    const char * lambdafile;
    const char * filterfile;
    if (!PyArg_ParseTuple(args, "s|sss",
                          &cfile,
                          &vfile, &lfile, &ffile))
        return NULL;
    /*sprintf(path,"%s/data/templates",getenv("KCORRECT_DIR"));*/
    printf("%s : v %s l %s f %s.\n",
           path, vfile, lfile, ffile);
    printf("-- coef %s !\n", cfile);

    Py_INCREF(Py_None);

    return Py_None;
    
}

PyDoc_STRVAR(kcorrect_fit_coeffs_doc,
" doc "
);

static PyMethodDef kcorrect_methods[] = {
    {"fit_coeffs", kcorrect_fit_coeffs, METH_VARARGS, kcorrect_fit_coeffs_doc},
    {NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module provides kcorrect"
);

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
PyInit_kcorrect(void)
{
	PyObject *m;
	m = PyModule_Create(&kcorrectmodule);
	if (m == NULL)
            return NULL;
        kcorrectError = PyErr_NewException("kcorrect.error", NULL, NULL);
        Py_INCREF(kcorrectError);
        PyModule_AddObject(m, "error", kcorrectError);
        return m;
}
