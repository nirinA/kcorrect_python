#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "structseq.h"
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <kcorrect.h>
#include "numpy/npy_3kcompat.h"
#include "numpy/arrayobject.h"

static PyObject *_ztransformError;


/* local prototype */
#define TOL 1.e-6

#if PY_VERSION_HEX < 0x03000000
#define PyStructSequence_GET_ITEM(op, i) PyTuple_GET_ITEM(op, i)
#endif

/* local function */
int omegas(float , float);
int omegas(float omega0, float omegal0)
{
    if ((fabs(omega0+omegal0-1.)<TOL) || (fabs(omegal0)<TOL && omega0<=1. && omega0>=0.)) {
        return (1);
    }
    else {
        return 0;
    }
}
static PyObject *
_ztransform_z2dm(PyObject *self, PyObject *args)
{
    float dm, z, omega0=0.3, omegal0=0.7;
    if (!PyArg_ParseTuple(args, "f|ff", &z, &omega0, &omegal0))
	return NULL;
    if (!omegas(omega0, omegal0)) {
        PyErr_SetString(_ztransformError, "cannot compute with these value of omega\n");
        return NULL;
    }
    dm = z2dm(z, omega0, omegal0);
    return Py_BuildValue("f", dm);
}

PyDoc_STRVAR(_ztransform_z2dm_doc,
"z2dm(z | omega0=0.3, omegal0=0.7) -> dm\n"
"Converts redshift z into distance modulus d.\n"
"INPUTS:\n"
"    z           redshift\n"
"OPTIONAL INPUTS:\n"
"    omega0      omega_matter to use (default: 0.3)\n"
"    omegal0     omega_lambda to use (default: 0.7)\n"
"OUTPUTS:\n"
"    dm          distance modulus");

static PyObject *
_ztransform_dm2z(PyObject *self, PyObject *args, PyObject *kwds)
{
    float dm, z, omega0=0.3, omegal0=0.7;
    static char *kwlist[] = {"dm", "omega0", "omegal0", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "f|ff", kwlist,
                                     &dm, &omega0, &omegal0))
        return NULL;
    if (!omegas(omega0, omegal0)) {
        PyErr_SetString(_ztransformError, "cannot compute with these value of omega\n");
        return NULL;
    }
    z = dm2z(dm, omega0, omegal0);
    return Py_BuildValue("f", z);
}

PyDoc_STRVAR(_ztransform_dm2z_doc,
"dm2z(dm | omega0=0.3, omegal0=0.7) -> z\n"
"Converts distance modulus into z\n"
"INPUTS:\n"
"    dm           distance modulus\n"
"OPTIONAL INPUTS:\n"
"    omega0      omega_matter to use (default: 0.3)\n"
"    omegal0     omega_lambda to use (default: 0.7)\n"
"OUTPUTS:\n"
"    z          redshift");

static PyObject *
_ztransform_ztoV(PyObject *self, PyObject *args, PyObject *kwds)
{
    float z, V, omega0=0.3, omegal0=0.7;
    static char *kwlist[] = {"z", "omega0", "omegal0", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "f|ff", kwlist,
                                     &z, &omega0, &omegal0))
        return NULL;
    if (!omegas(omega0, omegal0)) {
        PyErr_SetString(_ztransformError, "cannot compute with these value of omega\n");
        return NULL;
    }
    V = ztoV(z, omega0, omegal0);
    return Py_BuildValue("f", V);
}

PyDoc_STRVAR(_ztransform_ztoV_doc,
"ztoV(z | omega0=0.3, omegal0=0.7) -> V\n"
"INPUTS:\n"
"    z - redshift\n"
"OPTIONAL INPUTS:\n"
"    omega0 - matter density (default 0.3)\n"
"    omegal0 - vacuum energy density (default 0.7)\n"
"COMMENTS:\n"
"    returns R^3 in equation V=4*PI*R^3/3.");


static PyMethodDef _ztransform_methods[] = {
    {"z2dm", _ztransform_z2dm, METH_VARARGS, _ztransform_z2dm_doc},
    {"dm2z", (PyCFunction)_ztransform_dm2z, METH_VARARGS|METH_KEYWORDS, _ztransform_dm2z_doc},
    {"ztoV", (PyCFunction)_ztransform_ztoV, METH_VARARGS|METH_KEYWORDS, _ztransform_ztoV_doc},
    {NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module provides utilities\n\
to convert among redshift, lookback time,\n\
comoving distance, volume element, enclosed volume, distance\n\
modulus, and angular diameter distance.\n\
\n\
Works up to redshift 10 or so.\n\
\n\
adapted to Python from ztransform.c\n\
Original C code written by Blanton et al.");

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef _ztransformmodule = {
	PyModuleDef_HEAD_INIT,
	"_ztransform",
	module_doc,
	-1,
	_ztransform_methods,
	NULL,
	NULL,
	NULL,
	NULL
};

PyMODINIT_FUNC
PyInit__ztransform(void)
{
	PyObject *m;
	m = PyModule_Create(&_ztransformmodule);
	import_array();
	if (m == NULL)
            return NULL;
#else
PyMODINIT_FUNC
init_ztransform(void)
{
	PyObject *m;
	m = Py_InitModule3("_ztransform", _ztransform_methods, module_doc);
	import_array();
	if (m == NULL)
            goto finally;
#endif
        _ztransformError = PyErr_NewException("_ztransform.error", NULL, NULL);
        Py_INCREF(_ztransformError);
        PyModule_AddObject(m, "error", _ztransformError);

#if PY_VERSION_HEX >= 0x03000000
        return m;
#else
        finally:
        return;
#endif
}
