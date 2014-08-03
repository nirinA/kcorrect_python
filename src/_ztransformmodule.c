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
#if PY_VERSION_HEX < 0x03000000
#define PyStructSequence_GET_ITEM(op, i) PyTuple_GET_ITEM(op, i)
#endif

/* local function */
static PyObject *
_ztransform_z2dm(PyObject *self, PyObject *args, PyObject *kwds)
{
    float dm, z, omega0=0.3, omegal0=0.7;
    static char *kwlist[] = {"z", "omega0", "omegal0", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "f|ff", kwlist,
                                     &z, &omega0, &omegal0))
	return NULL;
    dm = z2dm(z, omega0, omegal0);
    return Py_BuildValue("f", dm);
}

PyDoc_STRVAR(_ztransform_z2dm_doc,
"(z, omega0=0.3, omegal0=0.7) -> dm\n\n"
"Converts redshift z into distance modulus dm.\n"
"INPUTS:\n"
"    z                redshifts\n"
"OPTIONAL INPUTS:\n"
"    omega0           omega_matter to use (default: 0.3)\n"
"    omegal0          omega_lambda to use (default: 0.7)\n"
"OUTPUTS:\n"
"    dm               distance modulus    ");


static PyMethodDef _ztransform_methods[] = {
    {"z2dm", (PyCFunction)_ztransform_z2dm, METH_VARARGS|METH_KEYWORDS, _ztransform_z2dm_doc},
    {NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(module_doc,
"This module provides utilities,\n\
utilities to convert among redshift, lookback time,\n\
comoving distance, volume element, enclosed volume, distance\n\
modulus, and angular diameter distance.\n\
\n\
Works up to redshift 10 or so.\n\
\n\
adapted to Python from ztransform.c\n\
Originally started circa 1999, Michael Blanton\n\
");

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
