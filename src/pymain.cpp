#if 1
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


typedef struct regMat {
	matrix<double> mat;
	matrix<double> refPts;
	matrix<double> movPts;
} regMat;

static PyObject *stackgreg_register(PyObject *self, PyObject *args);
static PyObject *stackgreg_transform(PyObject *self, PyObject *args);
//static PyObject *stackgreg_register_transform(PyObject *self, PyObject *args);
//static PyObject *stackgreg_get_points(PyObject *self, PyObject *args);
//static PyObject *stackgreg_get_matrix(PyObject *self, PyObject *args);

static char pystackreg_docs[] = "PyStackReg\n";
static PyMethodDef module_methods[] = {
    {"_register", (PyCFunction)stackgreg_register, METH_VARARGS, pystackreg_docs},
	{"_transform", (PyCFunction)stackgreg_transform, METH_VARARGS, pystackreg_docs},
	//{"register_transform", (PyCFunction)stackgreg_register_transform, METH_VARARGS, pystackreg_docs},
	//{"get_points", (PyCFunction)stackgreg_get_points, METH_VARARGS, pystackreg_docs},
	//{"get_matrix", (PyCFunction)stackgreg_get_matrix, METH_VARARGS, pystackreg_docs},
    {NULL}
};

static struct PyModuleDef pystackreg =
{
    PyModuleDef_HEAD_INIT,
    "pystackreg", /* name of module */
    "pystackreg\n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_stackreg(void)
{
    /* Load `numpy` functionality. */
    import_array();

    return PyModule_Create(&pystackreg);
}

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


	TurboRegTransform tform(movImg, movMsk, movPH, refImg, refMsk, refPH, transformation, false);

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


	TurboRegTransform tform(movImg, movMsk, movPH, transformation, false);

	std::vector<double> imgout = tform.doFinalTransform(movImg, m);
	return imgout;

}


PyObject *stackgreg_register(PyObject *self, PyObject *args) {

    PyObject *ref, *mov;
    regMat rm;
    unsigned char tf;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOB", &ref, &mov, &tf)) {
    	return NULL;
    }

    if ((tf != TRANSLATION) && (tf != RIGID_BODY) && (tf != SCALED_ROTATION) && (tf != AFFINE)) {
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

    registerImg(img_ref, img_mov, RIGID_BODY, Ny_ref, Nx_ref, rm); // width and height (Nx/Ny) have to be swapped!

    /* clean up */
    Py_XDECREF(ref_array);
    Py_XDECREF(mov_array);

    npy_intp dims_mat[2] = {rm.mat.nrows(), rm.mat.ncols()};
    npy_intp dims_pts[2] = {rm.refPts.nrows(), rm.refPts.ncols()};
	PyObject *retMat = PyArray_SimpleNew(2, (npy_intp*) &dims_mat, NPY_DOUBLE);
	PyObject *retPtsRef = PyArray_SimpleNew(2, (npy_intp*) &dims_pts, NPY_DOUBLE);
	PyObject *retPtsMov = PyArray_SimpleNew(2, (npy_intp*) &dims_pts, NPY_DOUBLE);

	memcpy((void*)PyArray_DATA(retMat),    rm.mat.begin(), (dims_mat[0] * dims_mat[1] * sizeof(double)));
	memcpy((void*)PyArray_DATA(retPtsRef), rm.refPts.begin(), (dims_pts[0] * dims_pts[1] * sizeof(double)));
	memcpy((void*)PyArray_DATA(retPtsMov), rm.movPts.begin(), (dims_pts[0] * dims_pts[1] * sizeof(double)));


    return Py_BuildValue("OOO", retMat, retPtsRef, retPtsMov);

}


PyObject *stackgreg_transform(PyObject *self, PyObject *args) {

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


//PyObject *stackgreg_register_transform(PyObject *self, PyObject *args) {return NULL;}
//PyObject *stackgreg_get_points(PyObject *self, PyObject *args) {return NULL;}
//PyObject *stackgreg_get_matrix(PyObject *self, PyObject *args) {return NULL;}

/*
static PyObject *std_std(PyObject *self, PyObject *args)
{
  PyObject* input;
  PyArg_ParseTuple(args, "O", &input);

  int size = PyList_Size(input);

  std::vector<double> list;
  list.resize(size);

  for(int i = 0; i < size; i++) {
    list[i] = PyFloat_AS_DOUBLE(PyList_GET_ITEM(input, i));
  }
  return PyFloat_FromDouble(standardDeviation(list));
}

static PyMethodDef std_methods[] = {
	{"standard_dev", std_standard_dev,	METH_VARARGS,
	 "Return the standard deviation of a list."},
	{NULL,		NULL}
};

extern void initstd(void)
{
	PyImport_AddModule("std");
	Py_InitModule("std", std_methods);
}

int main(int argc, char **argv)
{
	return 0;
}*/

#endif
