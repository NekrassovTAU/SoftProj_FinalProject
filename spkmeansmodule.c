#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "spkmeans.h"

#define ASSERT_ERROR(expr) \
    if (!(expr)){             \
        printf("%s", "An Error Has Occured"); \
        exit(0);\
    }

static PyObject *initializeCProcess(PyObject *self, PyObject *args);

static PyObject *KMeansPlusPlusIntegration(PyObject *self, PyObject *args);


static PyObject *initializeCProcess(PyObject *self, PyObject *args){
    /* argument in args is sys.argv from python only */

    /* Initialization and argument parsing */
    int i, j, argc, row_count, col_count;
    char **argv, *temp;
    double **ret_matrix;
    PyObject *argv_PyArray, *retPyMatrix, **tmpListPointer;
    if (!PyArg_ParseTuple(args, "O",
                          &argv_PyArray)){
        Py_DECREF(argv_PyArray);
        return NULL;
    }

    /* Parsing Python objects to C objects */
    argc = (int) PyList_Size(argv_PyArray);
    if(argc < 1){
        Py_DECREF(argv_PyArray);
        return NULL;
    }

    argv = calloc(argc, sizeof(char *));
    ASSERT_ERROR(argv != NULL)
    for (i = 0; i < argc; i++){
        temp = PyBytes_AsString(PyList_GetItem(argv_PyArray, i));
        argv[i] = strdup(temp);
    }

    ret_matrix = checkArgs(argc, argv, 1, &row_count, &col_count);

    /*Returning the data as a Python list of lists*/
    retPyMatrix = PyList_New(row_count);
    tmpListPointer = calloc(row_count, sizeof(PyObject*));
    for (i = 0; i < row_count ; i++){
        tmpListPointer[i] = PyList_New(col_count);
        for (j = 0; j < col_count ; j++){
            PyList_SetItem(tmpListPointer[i], j, PyFloat_FromDouble(ret_matrix[i][j]));
        }
        PyList_SetItem(retPyMatrix, i, tmpListPointer[i]);
    }

    /* Freeing allocated memory and finishing call*/
    for(i = 0; i < argc; i++){
        free(argv[i]);
    }
    free(argv);
    freeMatrix(&ret_matrix);

    for(i = 0; i < row_count ; i++){
        Py_DECREF(tmpListPointer[i]); //might be problematic? need to check
    }
    free(tmpListPointer);

    Py_DECREF(argv_PyArray);

    return retPyMatrix;
}

static PyObject *KMeansPlusPlusIntegration(PyObject *self, PyObject *args){
    /* argument in args is retPyMatrix from initializeCProcess,
     * and a list of initial centroids*/

    /* Initialization and argument parsing */
    int i, j, k, d, arraySize, *init_centroids;
    double **matrix, **retMatrix;
    PyObject *init_centroidsPyList, *matrixPyList, *retPyMatrix, *tmpPyList;

    if (!PyArg_ParseTuple(args, "OO",
                          &matrixPyList,
                          &init_centroidsPyList)){
        Py_DECREF(matrixPyList);
        Py_DECREF(init_centroidsPyList);
        return NULL;
    }

    /* Parsing Python objects to C objects */
    arraySize = (int) PyList_Size(matrixPyList);
    k = (int) PyList_Size(init_centroidsPyList);
    d = k;
    matrix = createMatrix(arraySize, d);
    init_centroids = calloc(k, sizeof(int));

    for (i = 0 ; i < arraySize ; i++){
        tmpPyList = PyList_GetItem(matrixPyList, i);
        for (j = 0 ; j < d; j++){
            matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(tmpPyList, j));
        }
    }

    for (i = 0; i < k; i++){
        init_centroids[i] = (int) PyLong_AsLong(PyList_GetItem(init_centroidsPyList, i));
    }

    /* Running KMeans algorithm in C */
    retMatrix = KMeansAlgorithm(k, &d, arraySize, &matrix, 1, &init_centroids);

    /* Parsing C Objects back to Python objects */
    retPyMatrix = PyList_New(arraySize * d);

    /* Freeing allocated memory and finishing call*/

    return NULL;
}

/* Defining the methods that are visible to the C extension spkmeansmodule*/
static PyMethodDef spkmeans_CAPI_Methods[] = {
        {"initializeCProcess", (PyCFunction) initializeCProcess, METH_VARARGS,
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