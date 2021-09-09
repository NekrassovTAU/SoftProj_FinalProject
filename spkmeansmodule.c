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

static PyObject ** createPyListArray(int n, int m);


static PyObject *
matrixToPyMatrix(double ***matrix, int row_count, int col_count);

/**
 * Receives arguments from Python, and processes them through the C-based
 * algorithm
 * Returns matrix containing results for printing / K-Means++
 * */
static PyObject *initializeCProcess(PyObject *self, PyObject *args){
    /* argument in args is sys.argv from python only */

    /* Initialization and argument parsing */
    int i, argc, row_count, col_count;
    char **argv, *temp;
    double **ret_matrix;
    PyObject *argv_PyArray, *retPyMatrix;

    if (!PyArg_ParseTuple(args, "O!", &PyList_Type,
                          &argv_PyArray)){
        Py_DECREF(argv_PyArray);
        Py_RETURN_NONE;
    }

    /* Parsing Python objects to C objects */
    argc = (int) PyObject_Length(argv_PyArray);
    if(argc < 0){
        Py_DECREF(argv_PyArray);
        printf("%s", "Invalid Input!\n");
        Py_RETURN_NONE;
    }

    /** TODO: Put in seperate function? */
    argv = calloc(argc, sizeof(char *));
    ASSERT_ERROR(argv != NULL)
    for (i = 0; i < argc; i++){
        temp = PyBytes_AsString(PyUnicode_AsEncodedString(
                PyList_GetItem(argv_PyArray, i),
                "utf-8", "strict"));
        argv[i] = strdup(temp);
    }

    ret_matrix = startCSPKMeansExec(argc, argv, 1, &row_count, &col_count);

    /*Returning the data as a Python list of lists*/
    retPyMatrix = matrixToPyMatrix(&ret_matrix, row_count, col_count);

    /* Freeing allocated memory and finishing call*/
    for(i = 0; i < argc; i++){
        free(argv[i]);
    }
    free(argv);
    freeMatrix(&ret_matrix);

    Py_DECREF(argv_PyArray);

    return retPyMatrix;
}

/**
 * Receives the matrix T and the result of K-Means++ (initial centroids),
 * runs it through the KMeans algorithm, and returns the final centroids
 * */
static PyObject *KMeansPlusPlusIntegration(PyObject *self, PyObject *args){
    /* argument in args is retPyMatrix from initializeCProcess,
     * and a list of initial centroids*/

    /* Initialization and argument parsing */
    int i, j, k, arraySize, *init_centroids;
    double **matrix, **retMatrix;
    PyObject *init_centroidsPyList, *matrixPyList, *retPyMatrix, *tmpPyList;

    if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type,
                          &matrixPyList, &PyList_Type,
                          &init_centroidsPyList)){
        Py_DECREF(matrixPyList);
        Py_DECREF(init_centroidsPyList);
        Py_RETURN_NONE;
    }

    /* Parsing Python objects to C objects */
    arraySize = (int) PyObject_Length(matrixPyList);
    k = (int) PyObject_Length(init_centroidsPyList);
    matrix = createMatrix(arraySize, k);
    init_centroids = calloc(k, sizeof(int));
    ASSERT_ERROR(init_centroids != NULL)

    tmpPyList = Py_None; /* initializing in order to not get warning from
                            * nova, works without it too*/

    for (i = 0 ; i < arraySize ; i++){
        tmpPyList = PyList_GetItem(matrixPyList, i);
        for (j = 0 ; j < k; j++){
            matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(tmpPyList, j));
        }
    }

    for (i = 0; i < k; i++){
        init_centroids[i] = (int)
                PyLong_AsLong(PyList_GetItem(init_centroidsPyList, i));
    }

    /* Running KMeans algorithm in C */
    retMatrix = KMeansAlgorithm(k, arraySize, &matrix, &init_centroids);

    /* Parsing C Objects back to Python objects */
    retPyMatrix = matrixToPyMatrix(&retMatrix, k, k);

    /* Freeing allocated memory and finishing call*/
    freeMatrix(&retMatrix);
    freeMatrix(&matrix);
    free(init_centroids);
    Py_DECREF(init_centroidsPyList);
    Py_DECREF(matrixPyList);
    Py_DECREF(tmpPyList);

    return retPyMatrix;
}


/**
 * Creates an array containing PyLists, used to initialize a python nested list
 * */
static PyObject ** createPyListArray(int n, int m){
    PyObject **array;
    int i;

    array = calloc(n, sizeof (PyObject* ));
    ASSERT_ERROR(array != NULL)

    for (i = 0; i < n; i++) {
        array[i] = PyList_New(m);
    }

    return array;
}

/**
 * Takes a C-matrix and converts it to a Python matrix (nested list, acts like
 * matrix)
 * @param matrix - C-matrix to be copied to python
 * @param row_count - amount of rows in matrix
 * @param col_count - amount of columns in matrix
 * @return
 */
static PyObject *
matrixToPyMatrix(double ***matrix, int row_count, int col_count) {

    int i, j;
    PyObject *retPyMatrix, **tmpList;

    retPyMatrix = PyList_New(row_count);
    tmpList = createPyListArray(row_count, col_count);
    for (i = 0; i < row_count ; i++){
        for (j = 0; j < col_count ; j++){
            PyList_SetItem(tmpList[i], j, PyFloat_FromDouble((*matrix)[i][j]));
        }
        PyList_SetItem(retPyMatrix, i, tmpList[i]);
    }

    return retPyMatrix;
}

/* Defining the methods that are visible to the C extension spkmeansmodule*/
static PyMethodDef spkmeans_CAPI_Methods[] = {
        {"initializeCProcess", (PyCFunction) initializeCProcess,
         METH_VARARGS,
         PyDoc_STR("SPKMeans Function implementation in C")},
        {"KMeansPlusPlusIntegration",
         (PyCFunction) KMeansPlusPlusIntegration,
         METH_VARARGS,
         PyDoc_STR("Uses KMeans++ results to determine spectral clusters")},
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