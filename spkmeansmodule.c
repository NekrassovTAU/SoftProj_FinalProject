#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>

#define ASSERT_ERROR(expr) \
    if (!(expr)){             \
        printf("%s", "An Error Has Occured"); \
        exit(0);\
    }

static PyObject *initProcess(PyObject *self, PyObject *args);

static PyObject *KMeansPlusPlusIntegration(PyObject *self, PyObject *args);


static PyObject *initProcess(PyObject *self, PyObject *args){
    //TODO: IMPLEMENT
    return NULL;
}

static PyObject *KMeansPlusPlusIntegration(PyObject *self, PyObject *args){
    //TODO: IMPLEMENT
    return NULL;
}

/* Defining the methods that are visible to the C extension spkmeansmodule*/
static PyMethodDef spkmeans_CAPI_Methods[] = {
        {"initProcess", (PyCFunction) initProcess, METH_VARARGS,
         PyDoc_STR("SPKMeans Function implementation in C")},
        {"KMeansPlusPlusIntegration", (PyCFunction) KMeansPlusPlusIntegration,
         METH_VARARGS,
         PyDoc_STR("Uses KMeans++ python results to determine spectral clusters")},
        {NULL, NULL, 0, NULL}
};

/* Defining the module spkmeansmodule*/
static struct PyModuleDef spkmeans_moduleDef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        spkmeans_CAPI_Methods,
        NULL,
        NULL,
        NULL,
        NULL
};

/* Building the module spkmeansmodule using definitions in spkmeans_moduleDef*/
PyMODINIT_FUNC
PyInit_spkmeansmodule (void)
{
    PyObject *m;
    m = PyModule_Create(&spkmeans_moduleDef);
    if(!m){
        return NULL;
    }
    return m;
}