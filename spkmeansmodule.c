#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include "string.h"
#include "stdlib.h"
#include "spkmeans.h"
#include "math.h"

/*decleraton of global variables*/
int k;
int iter;
double epsilon;
int numPoints;
int pointLength;
double **points;
double **means;
double **finalMeans;
double **dataPoints;

/*prototypes*/
static PyObject* spk(PyObject *self, PyObject *args);
static PyObject* wam(PyObject* self, PyObject* args);
static PyObject* ddg(PyObject* self, PyObject* args);
static PyObject* gl(PyObject* self, PyObject* args);
static PyObject* jacobi(PyObject* self, PyObject* args);
static PyObject* pre_spk(PyObject* self, PyObject* args);
PyObject *convertToPyObject(int dim_row, int dim_col, double **Matrix);


/*return final centoids*/
static PyObject* spk(PyObject* self, PyObject* args){
    PyObject* py_points;
    PyObject* py_means;
    int i;

     /*parses the Python arguments into C arguments*/
    if(!PyArg_ParseTuple(args, "iidOiiO", &k, &iter, &epsilon, &py_points, &numPoints, &pointLength, &py_means )) {
        return NULL; 
    }
    /*update points array from the pyObject*/
    points=initial2Darray(numPoints, pointLength);
    for (i = 0; i < numPoints; i++) {
        int j;
        for (j = 0; j < pointLength; j++) {
            points[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_points, i), j));
        } 
       
    }
    /*update means array from the pyObject*/
    means=initial2Darray(k, pointLength);
    for (i = 0; i < k; i++){
        int l;
        for (l = 0; l < pointLength; l++) {
            means[i][l] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_means, i), l));
        }
    }

    double** finalMeans=k_means(k, iter, epsilon, numPoints, pointLength, means, points);
    PyObject* p= convertToPyObject(k, pointLength, finalMeans); /*convert finalMeans into pyObject*/
    freeMemory_kmeans();
    return p;

}

/*returns W*/
static PyObject* wam(PyObject* self, PyObject* args){
    int i;
    PyObject* py_points;

    /*parses the Python arguments into C arguments*/
    if(!PyArg_ParseTuple(args, "Oii", &py_points, &numPoints, &pointLength )) {
        return NULL; 
    }

    /*allocating space for dataPoints which stores the points and insert values to it*/
    dataPoints = calloc(numPoints,sizeof(double *));
    if (dataPoints == NULL){
        freeMemoryAndExit();
    }
    for (i = 0; i < numPoints; i++) {
        dataPoints[i] = malloc(pointLength* sizeof(double));
        if (dataPoints[i] == NULL) {
             freeMemoryAndExit();
        }
        int j;
        for (j = 0; j < pointLength; j++) {
            dataPoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_points, i), j));
        } 
    }

    double** W=wam_calc(dataPoints, numPoints, pointLength);
    PyObject* p= convertToPyObject(numPoints, numPoints, W); /*convert W into pyObject*/
    freeMemory();
    return p;
}

/*returns D*/
static PyObject* ddg(PyObject* self, PyObject* args){
    int i;
    PyObject* py_points;

    /*parses the Python arguments into C arguments*/
    if(!PyArg_ParseTuple(args, "Oii", &py_points, &numPoints, &pointLength )) {
        return NULL; 
    }

    /*allocating space for dataPoints which stores the points and insert values to it*/
    dataPoints = calloc(numPoints,sizeof(double *));
    if (dataPoints == NULL){
        freeMemoryAndExit();
    }
    for (i = 0; i < numPoints; i++) {
        dataPoints[i] = malloc(pointLength* sizeof(double));
        if (dataPoints[i] == NULL) {
             freeMemoryAndExit();
        }
        int j;
        for (j = 0; j < pointLength; j++) {
            dataPoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_points, i), j));
        } 
    }

    double** D=ddg_calc(dataPoints, numPoints, pointLength);
    PyObject* p= convertToPyObject(numPoints, numPoints, D); /*convert D into pyObject*/
    freeMemory();
    return p;
}

static PyObject* gl(PyObject* self, PyObject* args){
    int i;
    PyObject* py_points;

    /*parses the Python arguments into C arguments*/
    if(!PyArg_ParseTuple(args, "Oii", &py_points, &numPoints, &pointLength )) {
        return NULL; 
    }

    /*allocating space for dataPoints which stores the points and insert values to it*/
    dataPoints = calloc(numPoints,sizeof(double *));
    if (dataPoints == NULL){
        freeMemoryAndExit();
    }
    for (i = 0; i < numPoints; i++) {
        dataPoints[i] = malloc(pointLength* sizeof(double));
        if (dataPoints[i] == NULL) {
             freeMemoryAndExit();
        }
        int j;
        for (j = 0; j < pointLength; j++) {
            dataPoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_points, i), j));
        } 
    }

    double** L=gl_calc(dataPoints, numPoints, pointLength);
    PyObject* p = convertToPyObject(numPoints, numPoints, L); /*convert L into pyObject*/
    freeMemory();
    return p;
}

static PyObject* jacobi(PyObject* self, PyObject* args){
    int i;
    PyObject* py_points;

    /*parses the Python arguments into C arguments*/
    if(!PyArg_ParseTuple(args, "Oii", &py_points, &numPoints, &pointLength )) {
        return NULL; 
    }

    /*allocating space for dataPoints which stores the points and insert values to it*/
    dataPoints = calloc(numPoints,sizeof(double *));
    if (dataPoints == NULL){
        freeMemoryAndExit();
    }
    for (i = 0; i < numPoints; i++) {
        dataPoints[i] = malloc(pointLength* sizeof(double));
        if (dataPoints[i] == NULL) {
             freeMemoryAndExit();
        }
        int j;
        for (j = 0; j < pointLength; j++) {
            dataPoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_points, i), j));
        } 
    }
    double** JACOBI=jacobi_calc(dataPoints, numPoints, pointLength);
    PyObject* p = convertToPyObject(numPoints+1, numPoints, JACOBI); /*convert JACOBI into pyObject*/
    freeMemory();
    return p;
}

static PyObject* pre_spk(PyObject* self, PyObject* args){
    int i;
    int arg_k;
    PyObject* py_points;

    /*parses the Python arguments into C arguments*/
    if(!PyArg_ParseTuple(args, "Oiii", &py_points, &arg_k, &numPoints, &pointLength )) {
        return NULL; 
    }

    /*allocating space for dataPoints which stores the points and insert values to it*/
    dataPoints = calloc(numPoints,sizeof(double *));
    if (dataPoints == NULL){
        freeMemoryAndExit();
    }
    for (i = 0; i < numPoints; i++) {
        dataPoints[i] = malloc(pointLength* sizeof(double));
        if (dataPoints[i] == NULL) {
             freeMemoryAndExit();
        }
        int j;
        for (j = 0; j < pointLength; j++) {
            dataPoints[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(py_points, i), j));
        } 
    }

    double** SPK=spk_calc(dataPoints, arg_k, numPoints, pointLength);
    k= SPK[numPoints][0]; /*get k value from SPK matrix*/
    PyObject* p= convertToPyObject(numPoints, k, SPK); /*convert SPK (without the last row) into pyObject*/
    freeMemory();
    return p;
}

static PyMethodDef kmeansMethods[] = {
        {"spk", (PyCFunction) spk, METH_VARARGS,  PyDoc_STR("calculates K-Means")}, 
        {"pre_spk", (PyCFunction) pre_spk, METH_VARARGS,  PyDoc_STR("calculates first k eigenvectors")},
        {"wam", (PyCFunction) wam, METH_VARARGS,  PyDoc_STR("calculates W")},
        {"ddg", (PyCFunction) ddg, METH_VARARGS,  PyDoc_STR("calculates D")},
        {"gl", (PyCFunction) gl, METH_VARARGS,  PyDoc_STR("calculates L")},
        {"jacobi", (PyCFunction) jacobi, METH_VARARGS,  PyDoc_STR("calculates eigenvalues and eigenvectors")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mykmeanssp = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp", /* name of module */
        NULL, /* module documentation, may be NULL */
        -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        kmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}

/*convers C matrix to pyObject*/
PyObject *convertToPyObject(int dim_row, int dim_col, double **Matrix){
    int i, j;
    PyObject *py_list = NULL;
    PyObject *py_sublist = NULL;
    py_list = PyList_New(0);
    for (i = 0; i < dim_row; i++){
        py_sublist = PyList_New(0);
        for (j = 0; j < dim_col; j++){
            PyList_Append(py_sublist, PyFloat_FromDouble(Matrix[i][j]));
        }
        PyList_Append(py_list, py_sublist);
    }
    return py_list;
}