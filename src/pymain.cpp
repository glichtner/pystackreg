#include <Python.h>
#include <numpy/arrayobject.h>

#include <vector>
#include <numeric>
#include <iterator>


static PyObject *stackreg_std(PyObject *self, PyObject *args);
static PyObject *stackreg_mean(PyObject *self, PyObject *args);
static PyObject *stackreg_mean_center(PyObject *self, PyObject *args);

static char pystackreg_docs[] = "PyStackReg\n";
static PyMethodDef module_methods[] = {
    {"std", (PyCFunction)stackreg_std, METH_VARARGS, pystackreg_docs},
	{"mean", (PyCFunction)stackreg_mean, METH_VARARGS, pystackreg_docs},
	{"mean_center", (PyCFunction)stackreg_mean_center, METH_VARARGS, pystackreg_docs},
    {NULL}
};

static struct PyModuleDef PyStackReg =
{
    PyModuleDef_HEAD_INIT,
    "PyStackReg", /* name of module */
    "PyStackReg\n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_PyStackReg(void)
{
    /* Load `numpy` functionality. */
    import_array();

    return PyModule_Create(&PyStackReg);
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




PyObject *stackreg_mean(PyObject *self, PyObject *args) {
    PyObject *x_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "O", &x_obj))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);


    /* If that didn't work, throw an exception. */
    if (x_array == NULL) {
        Py_XDECREF(x_array);
        return NULL;
    }

    /* How many data points are there? */
    int N = (int)PyArray_DIM(x_array, 0);

    /* Get pointers to the data as C-types. */
    double *x    = (double*)PyArray_DATA(x_array);

	return PyFloat_FromDouble(mean(x, N));
}

PyObject *stackreg_mean_center(PyObject *self, PyObject *args) {

    PyObject *x_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "O", &x_obj))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);


    /* If that didn't work, throw an exception. */
    if (x_array == NULL) {
        Py_XDECREF(x_array);
        return NULL;
    }

    /* How many data points are there? */
    int N = (int)PyArray_DIM(x_array, 0);
    npy_intp dims[1] = {N};

    /* Get pointers to the data as C-types. */
    double *x    = (double*)PyArray_DATA(x_array);


	PyObject *ret = PyArray_SimpleNew(1, (npy_intp*) &dims, NPY_DOUBLE);

	double avg = mean(x, N);

	memcpy((void*)PyArray_DATA(ret), x, (N * sizeof(double)));

	double* y = (double*)PyArray_DATA(ret);

	for(int i=0; i < N; i++) {
		y[i] -= avg;
	}


	return ret;
}

static PyObject *stackreg_std(PyObject *self, PyObject *args)
{
    PyObject *x_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "O", &x_obj))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);


    /* If that didn't work, throw an exception. */
    if (x_array == NULL) {
        Py_XDECREF(x_array);
        return NULL;
    }

    /* How many data points are there? */
    int N = (int)PyArray_DIM(x_array, 0);

    /* Get pointers to the data as C-types. */
    double *x    = (double*)PyArray_DATA(x_array);

    /* Call the external C function to compute the chi-squared. */
    double value = cstd(x, N);

    /* Clean up. */
    Py_DECREF(x_array);

    if (value < 0.0) {
        PyErr_SetString(PyExc_RuntimeError,
                    "Chi-squared returned an impossible value.");
        return NULL;
    }

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}


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
