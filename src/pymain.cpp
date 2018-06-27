#if 0
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
#include "matrix.h"

static PyObject *stackgreg_register(PyObject *self, PyObject *args);
static PyObject *stackgreg_transform(PyObject *self, PyObject *args);
static PyObject *stackgreg_register_transform(PyObject *self, PyObject *args);
static PyObject *stackgreg_get_points(PyObject *self, PyObject *args);
static PyObject *stackgreg_get_matrix(PyObject *self, PyObject *args);

static char pystackreg_docs[] = "PyStackReg\n";
static PyMethodDef module_methods[] = {
    {"register", (PyCFunction)stackgreg_register, METH_VARARGS, pystackreg_docs},
	{"transform", (PyCFunction)stackgreg_transform, METH_VARARGS, pystackreg_docs},
	{"register_transform", (PyCFunction)stackgreg_register_transform, METH_VARARGS, pystackreg_docs},
	{"get_points", (PyCFunction)stackgreg_get_points, METH_VARARGS, pystackreg_docs},
	{"get_matrix", (PyCFunction)stackgreg_get_matrix, METH_VARARGS, pystackreg_docs},
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


double mean(double *x, int N) {
	double avg = 0;

	for(int i=0; i < N; i++) {
		avg += x[i];
	}

	return avg / (double)N;
}

double cstd(double *x, int N) {

	double avg = mean(x, N);
	double std = 0;

	for(int i=0; i < N; i++) {
		std += std::pow(x[i] - avg, 2);
	}

	return std / ((double)N - 1);
}

std::vector<double> registerImg(double *pDataRef, double *pDataMov, Transformation transformation, int width, int height) {
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

	std::vector<double> imgout = tform.doFinalTransform(width, height);
	return imgout;
}


PyObject *stackgreg_register(PyObject *self, PyObject *args) {

    PyObject *ref, *mov;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &ref, &mov)) {
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


    npy_intp dims[2] = {Nx_ref, Ny_ref};

    /* Get pointers to the data as C-types. */
    double *img_ref    = (double*)PyArray_DATA(ref_array);
    double *img_mov    = (double*)PyArray_DATA(mov_array);

    std::vector<double> imgout = registerImg(img_ref, img_mov, RIGID_BODY, Ny_ref, Nx_ref); // width and height (Nx/Ny) have to be swapped!

	PyObject *ret = PyArray_SimpleNew(2, (npy_intp*) &dims, NPY_DOUBLE);

	memcpy((void*)PyArray_DATA(ret), &imgout[0], (Nx_ref * Ny_ref * sizeof(double)));

    /* clean up */
    Py_XDECREF(ref_array);
    Py_XDECREF(mov_array);


	return ret;
}
PyObject *stackgreg_transform(PyObject *self, PyObject *args) {return NULL;}
PyObject *stackgreg_register_transform(PyObject *self, PyObject *args) {return NULL;}
PyObject *stackgreg_get_points(PyObject *self, PyObject *args) {return NULL;}
PyObject *stackgreg_get_matrix(PyObject *self, PyObject *args) {return NULL;}

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
