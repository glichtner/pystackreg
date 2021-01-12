/*
 * pymain.cpp
 *
 * C++/Python Interface functions for the PyStackReg package.
 *
 *  Created on: Jun 27, 2018
 *      Author: Gregor Lichtner
 */

#include <Python.h>
#include <numpy/arrayobject.h>

#include <vector>
#include <numeric>
#include <iterator>

#include "TurboReg.h"
#include "TurboRegImage.h"
#include "TurboRegMask.h"
#include "TurboRegPointHandler.h"
#include "TurboRegTransform.h"

#if PY_MAJOR_VERSION >= 3
#define ISPY3
#endif


typedef struct regMat {
	matrix<double> mat;
	matrix<double> refPts;
	matrix<double> movPts;
} regMat;

static PyObject *turbogreg_register(PyObject *self, PyObject *args);
static PyObject *turbogreg_transform(PyObject *self, PyObject *args);



static char pystackreg_docs[] = "PyStackReg\n";
static PyMethodDef module_methods[] = {
    {"_register", (PyCFunction)turbogreg_register, METH_VARARGS, pystackreg_docs},
	{"_transform", (PyCFunction)turbogreg_transform, METH_VARARGS, pystackreg_docs},
    {NULL}
};

#ifdef ISPY3
static struct PyModuleDef pystackreg =
{
    PyModuleDef_HEAD_INIT,
    "turboreg", /* name of module */
    "turboreg\n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};


PyMODINIT_FUNC PyInit_turboreg(void)
{
    /* Load `numpy` functionality. */
    import_array();

    return PyModule_Create(&pystackreg);
}
#else

PyMODINIT_FUNC initturboreg(void)
{
    /* Load `numpy` functionality. */
    import_array();
    Py_InitModule("turboreg", module_methods);
}

#endif

bool registerImg(double *pDataRef, double *pDataMov, Transformation transformation, int width, int height, regMat &rm) {
	TurboRegImage refImg(pDataRef, width, height, transformation, true);
	TurboRegImage movImg(pDataMov, width, height, transformation, false);

	TurboRegPointHandler refPH(refImg, transformation);
	TurboRegPointHandler movPH(movImg, transformation);

	TurboRegMask refMsk(refImg);
	TurboRegMask movMsk(movImg);

	refMsk.clearMask();
	movMsk.clearMask();

	int pyramidDepth = getPyramidDepth(
			movImg.getWidth(), movImg.getHeight(),
			refImg.getWidth(), refImg.getHeight()
			);
	refImg.setPyramidDepth(pyramidDepth);
	refMsk.setPyramidDepth(pyramidDepth);
	movImg.setPyramidDepth(pyramidDepth);
	movMsk.setPyramidDepth(pyramidDepth);

	refImg.init();
	refMsk.init();
	movImg.init();
	movMsk.init();


	TurboRegTransform tform(&movImg, &movMsk, &movPH, &refImg, &refMsk, &refPH, transformation, false);

	tform.doRegistration();


	rm.mat = tform.getTransformationMatrix();
	rm.refPts = refPH.getPoints();
	rm.movPts = movPH.getPoints();

	return true;
}

std::vector<double> transformImg(matrix<double> m, double *pDataMov, int width, int height) {

	Transformation transformation = getTransformationFromMatrix(m);

	TurboRegImage movImg(pDataMov, width, height, transformation, false);

	TurboRegPointHandler movPH(movImg, transformation);

	TurboRegMask movMsk(movImg);

	movMsk.clearMask();

	int pyramidDepth = getPyramidDepth(
			movImg.getWidth(), movImg.getHeight(),
			movImg.getWidth(), movImg.getHeight()
			);
	movImg.setPyramidDepth(pyramidDepth);
	movMsk.setPyramidDepth(pyramidDepth);

	movImg.init();
	movMsk.init();


	TurboRegTransform tform(&movImg, &movMsk, &movPH, transformation, false);

	std::vector<double> imgout = tform.doFinalTransform(&movImg, m);
	return imgout;

}


PyObject *turbogreg_register(PyObject *self, PyObject *args) {

    PyObject *ref, *mov;
    regMat rm;
    unsigned char tf;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOB", &ref, &mov, &tf)) {
    	return NULL;
    }

    if ((tf != TRANSLATION) && (tf != RIGID_BODY) && (tf != SCALED_ROTATION) && (tf != AFFINE) && (tf != BILINEAR)) {
    	PyErr_SetString(PyExc_ValueError, "Invalid transformation");
        return NULL;
    }

    /* Interpret the input objects as numpy arrays. */
    PyObject *ref_array = PyArray_FROM_OTF(ref, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mov_array = PyArray_FROM_OTF(mov, NPY_DOUBLE, NPY_IN_ARRAY);


    /* If that didn't work, throw an exception. */
    if (ref_array == NULL || mov_array == NULL) {
        Py_XDECREF(ref_array);
        Py_XDECREF(mov_array);
        return NULL;
    }

    int ndim_ref = (int)PyArray_NDIM(ref_array);
    int ndim_mov = (int)PyArray_NDIM(mov_array);

    if (ndim_ref != 2 || ndim_mov != 2) {
        Py_XDECREF(ref_array);
        Py_XDECREF(mov_array);
    	PyErr_SetString(PyExc_ValueError, "Input arrays must be two dimensional");
        return NULL;
    }

    /* How many data points are there? */
    int Nx_ref = (int)PyArray_DIM(ref_array, 0);
    int Ny_ref = (int)PyArray_DIM(ref_array, 1);
    int Nx_mov = (int)PyArray_DIM(mov_array, 0);
    int Ny_mov = (int)PyArray_DIM(mov_array, 1);

    if( Nx_ref != Nx_mov || Ny_ref != Ny_mov) {
        Py_XDECREF(ref_array);
        Py_XDECREF(mov_array);
    	PyErr_SetString(PyExc_ValueError, "Input arrays must of the same shape");
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *img_ref    = (double*)PyArray_DATA(ref_array);
    double *img_mov    = (double*)PyArray_DATA(mov_array);

    registerImg(img_ref, img_mov, (Transformation)tf, Ny_ref, Nx_ref, rm); // width and height (Nx/Ny) have to be swapped!

    /* clean up */
    Py_XDECREF(ref_array);
    Py_XDECREF(mov_array);

    npy_intp dims_mat[2];
    npy_intp dims_pts[2];

    dims_mat[0] = rm.mat.nrows();
    dims_mat[1] = rm.mat.ncols();
    dims_pts[0] = rm.refPts.nrows();
    dims_pts[1] = rm.refPts.ncols();

	PyObject *retMat = PyArray_SimpleNew(2, (npy_intp*) &dims_mat, NPY_DOUBLE);
	PyObject *retPtsRef = PyArray_SimpleNew(2, (npy_intp*) &dims_pts, NPY_DOUBLE);
	PyObject *retPtsMov = PyArray_SimpleNew(2, (npy_intp*) &dims_pts, NPY_DOUBLE);

	memcpy((void*)PyArray_DATA(retMat),    rm.mat.begin(), (dims_mat[0] * dims_mat[1] * sizeof(double)));
	memcpy((void*)PyArray_DATA(retPtsRef), rm.refPts.begin(), (dims_pts[0] * dims_pts[1] * sizeof(double)));
	memcpy((void*)PyArray_DATA(retPtsMov), rm.movPts.begin(), (dims_pts[0] * dims_pts[1] * sizeof(double)));


    return Py_BuildValue("OOO", retMat, retPtsRef, retPtsMov);

}


PyObject *turbogreg_transform(PyObject *self, PyObject *args) {

    PyObject *mov;
    PyObject *mat;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &mov, &mat)) {
    	return NULL;
    }

    /* Interpret the input objects as numpy arrays. */
    PyObject *mov_array = PyArray_FROM_OTF(mov, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mat_array = PyArray_FROM_OTF(mat, NPY_DOUBLE, NPY_IN_ARRAY);


    /* If that didn't work, throw an exception. */
    if (mov_array == NULL || mat_array == NULL) {
        Py_XDECREF(mat_array);
        Py_XDECREF(mov_array);
        return NULL;
    }

    int ndim_mov = (int)PyArray_NDIM(mov_array);
    int ndim_mat = (int)PyArray_NDIM(mat_array);

    if (ndim_mov != 2 || ndim_mat != 2) {
        Py_XDECREF(mat_array);
        Py_XDECREF(mov_array);
    	PyErr_SetString(PyExc_ValueError, "Input arrays must be two dimensional");
        return NULL;
    }

    /* How many data points are there? */
    int Nx_mov = (int)PyArray_DIM(mov_array, 0);
    int Ny_mov = (int)PyArray_DIM(mov_array, 1);
    int Nx_mat = (int)PyArray_DIM(mat_array, 0);
    int Ny_mat = (int)PyArray_DIM(mat_array, 1);


    if( Nx_mat != 2 || (
    		(Ny_mat != 1) && (Ny_mat != 3) && (Ny_mat != 4))) {
    	Py_XDECREF(mov_array);
    	Py_XDECREF(mat_array);
    	PyErr_SetString(PyExc_ValueError, "Transformation matrix must be of shape (2,1), (2,3) or (2,4)");
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *img_mov    = (double*)PyArray_DATA(mov_array);
    double *tmat    = (double*)PyArray_DATA(mat_array);

    matrix<double> m(Nx_mat, Ny_mat);
    memcpy(m.begin(), tmat, (Nx_mat * Ny_mat * sizeof(double)));



    std::vector<double> imgout = transformImg(m, img_mov, Ny_mov, Nx_mov); // width and height (Nx/Ny) have to be swapped!

    /* clean up */
    Py_XDECREF(mat_array);
    Py_XDECREF(mov_array);

	npy_intp dims[2] = {Nx_mov, Ny_mov};

	PyObject *ret = PyArray_SimpleNew(2, (npy_intp*) &dims, NPY_DOUBLE);

	memcpy((void*)PyArray_DATA(ret), &imgout[0], (Nx_mov * Ny_mov * sizeof(double)));


    return ret;


}
