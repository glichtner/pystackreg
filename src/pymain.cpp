/*
 * pymain.cpp
 *
 * C++/Python Interface functions for the PyStackReg package.
 *
 *  Created on: Jun 27, 2018
 *      Author: Gregor Lichtner
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

#include <climits>
#include <cmath>
#include <exception>
#include <iterator>
#include <new>
#include <numeric>
#include <string>
#include <vector>

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

/* Outcome of a native computation that ran without the interpreter lock. */
enum NativeStatus { NATIVE_OK, NATIVE_NO_MEMORY, NATIVE_ERROR };

/*
 * Largest coordinate magnitude (in pixels) that a transformation may produce.
 * The resampling code converts coordinates to int; staying well inside the
 * int range keeps those conversions defined.
 */
static const double MAX_COORDINATE = 1.0e9;

/* Set a Python exception for a failed native computation. */
static void setNativeError(NativeStatus status, const std::string &message) {
  if (status == NATIVE_NO_MEMORY) {
    PyErr_NoMemory();
  } else {
    PyErr_SetString(PyExc_RuntimeError, message.empty()
                                            ? "Unknown error in TurboReg"
                                            : message.c_str());
  }
}

/* Check that an image shape can be handled by the int-indexed native code. */
static bool checkImageSize(npy_intp rows, npy_intp cols, npy_intp minSize,
                           const char *what) {
  if (rows < minSize || cols < minSize) {
    PyErr_Format(PyExc_ValueError,
                 "%s must be at least %d x %d pixels, got %zd x %zd", what,
                 (int)minSize, (int)minSize, (Py_ssize_t)rows,
                 (Py_ssize_t)cols);
    return false;
  }
  if (rows > INT_MAX / 2 || cols > INT_MAX / 2 || rows > INT_MAX / cols) {
    PyErr_Format(PyExc_ValueError,
                 "%s is too large (%zd x %zd pixels); the number of pixels "
                 "must not exceed %d",
                 what, (Py_ssize_t)rows, (Py_ssize_t)cols, INT_MAX);
    return false;
  }
  return true;
}

static PyObject *turbogreg_register(PyObject *self, PyObject *args);
static PyArrayObject *turbogreg_transform(PyObject *self, PyObject *args);

static char pystackreg_docs[] = "PyStackReg\n";
static PyMethodDef module_methods[] = {
    {"_register", (PyCFunction)turbogreg_register, METH_VARARGS,
     pystackreg_docs},
    {"_transform", (PyCFunction)turbogreg_transform, METH_VARARGS,
     pystackreg_docs},
    {NULL}};

#ifdef ISPY3
static struct PyModuleDef pystackreg = {
    PyModuleDef_HEAD_INIT, "turboreg", /* name of module */
    "turboreg\n",                      /* module documentation, may be NULL */
    -1, /* size of per-interpreter state of the module, or -1 if the module
           keeps state in global variables. */
    module_methods};

PyMODINIT_FUNC PyInit_turboreg(void) {
  /* Load `numpy` functionality. */
  import_array();

  PyObject *module = PyModule_Create(&pystackreg);
#ifdef Py_GIL_DISABLED
  /* The native code keeps all of its state in objects local to one call, so
     the module does not need the GIL on free-threaded builds. */
  if (module != NULL) {
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
  }
#endif
  return module;
}
#else

PyMODINIT_FUNC initturboreg(void) {
  /* Load `numpy` functionality. */
  import_array();
  Py_InitModule("turboreg", module_methods);
}

#endif

bool registerImg(double *pDataRef, double *pDataMov,
                 Transformation transformation, int width, int height,
                 regMat &rm) {
  TurboRegImage refImg(pDataRef, width, height, transformation, true);
  TurboRegImage movImg(pDataMov, width, height, transformation, false);

  TurboRegPointHandler refPH(refImg, transformation);
  TurboRegPointHandler movPH(movImg, transformation);

  TurboRegMask refMsk(refImg);
  TurboRegMask movMsk(movImg);

  refMsk.clearMask();
  movMsk.clearMask();

  int pyramidDepth = getPyramidDepth(movImg.getWidth(), movImg.getHeight(),
                                     refImg.getWidth(), refImg.getHeight());
  refImg.setPyramidDepth(pyramidDepth);
  refMsk.setPyramidDepth(pyramidDepth);
  movImg.setPyramidDepth(pyramidDepth);
  movMsk.setPyramidDepth(pyramidDepth);

  refImg.init();
  refMsk.init();
  movImg.init();
  movMsk.init();

  TurboRegTransform tform(&movImg, &movMsk, &movPH, &refImg, &refMsk, &refPH,
                          transformation, false);

  tform.doRegistration();

  rm.mat = tform.getTransformationMatrix();
  rm.refPts = refPH.getPoints();
  rm.movPts = movPH.getPoints();

  return true;
}

std::vector<double> transformImg(matrix<double> m, double *pDataMov, int width,
                                 int height) {
  Transformation transformation = getTransformationFromMatrix(m);

  /* Resampling with a known transformation only needs the B-spline
     coefficients of the image; masks, gradients and pyramids are used
     exclusively by the registration and are not built here. */
  TurboRegImage movImg(pDataMov, width, height, transformation, false);
  movImg.initCoefficients();

  TurboRegTransform tform(&movImg, NULL, NULL, transformation, false);

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

  if ((tf != TRANSLATION) && (tf != RIGID_BODY) && (tf != SCALED_ROTATION) &&
      (tf != AFFINE) && (tf != BILINEAR)) {
    PyErr_SetString(PyExc_ValueError, "Invalid transformation");
    return NULL;
  }

  /* Interpret the input objects as numpy arrays. */
  PyObject *po_ref_array = PyArray_FROM_OTF(ref, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (po_ref_array == NULL) {
    return NULL;
  }
  PyObject *po_mov_array = PyArray_FROM_OTF(mov, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  PyArrayObject *ref_array = (PyArrayObject *)po_ref_array;
  PyArrayObject *mov_array = (PyArrayObject *)po_mov_array;

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

  /* the gradient filters read two samples per image line */
  if (!checkImageSize(PyArray_DIM(ref_array, 0), PyArray_DIM(ref_array, 1), 2,
                      "Images passed to the registration") ||
      !checkImageSize(PyArray_DIM(mov_array, 0), PyArray_DIM(mov_array, 1), 2,
                      "Images passed to the registration")) {
    Py_XDECREF(ref_array);
    Py_XDECREF(mov_array);
    return NULL;
  }

  /* How many data points are there? */
  int Nx_ref = (int)PyArray_DIM(ref_array, 0);
  int Ny_ref = (int)PyArray_DIM(ref_array, 1);
  int Nx_mov = (int)PyArray_DIM(mov_array, 0);
  int Ny_mov = (int)PyArray_DIM(mov_array, 1);

  if (Nx_ref != Nx_mov || Ny_ref != Ny_mov) {
    Py_XDECREF(ref_array);
    Py_XDECREF(mov_array);
    PyErr_SetString(PyExc_ValueError, "Input arrays must of the same shape");
    return NULL;
  }

  /* Get pointers to the data as C-types. */
  double *img_ref = (double *)PyArray_DATA(ref_array);
  double *img_mov = (double *)PyArray_DATA(mov_array);

  NativeStatus status = NATIVE_OK;
  std::string message;

  Py_BEGIN_ALLOW_THREADS
  try {
    registerImg(img_ref, img_mov, (Transformation)tf, Ny_ref, Nx_ref,
                rm); // width and height (Nx/Ny) have to be swapped!
  } catch (const std::bad_alloc &) {
    status = NATIVE_NO_MEMORY;
  } catch (const std::exception &e) {
    status = NATIVE_ERROR;
    try {
      message = e.what();
    } catch (...) {
    }
  } catch (...) {
    status = NATIVE_ERROR;
  }
  Py_END_ALLOW_THREADS

  /* clean up */
  Py_XDECREF(ref_array);
  Py_XDECREF(mov_array);

  if (status != NATIVE_OK) {
    setNativeError(status, message);
    return NULL;
  }

  npy_intp dims_mat[2];
  npy_intp dims_pts[2];

  dims_mat[0] = rm.mat.nrows();
  dims_mat[1] = rm.mat.ncols();
  dims_pts[0] = rm.refPts.nrows();
  dims_pts[1] = rm.refPts.ncols();

  PyArrayObject *retMat = (PyArrayObject*)PyArray_SimpleNew(2, (npy_intp *)&dims_mat, NPY_DOUBLE);
  PyArrayObject *retPtsRef = (PyArrayObject*)PyArray_SimpleNew(2, (npy_intp *)&dims_pts, NPY_DOUBLE);
  PyArrayObject *retPtsMov = (PyArrayObject*)PyArray_SimpleNew(2, (npy_intp *)&dims_pts, NPY_DOUBLE);

  if (retMat == NULL || retPtsRef == NULL || retPtsMov == NULL) {
    Py_XDECREF(retMat);
    Py_XDECREF(retPtsRef);
    Py_XDECREF(retPtsMov);
    return NULL;
  }

  memcpy((void *)PyArray_DATA(retMat), rm.mat.begin(),
         (dims_mat[0] * dims_mat[1] * sizeof(double)));
  memcpy((void *)PyArray_DATA(retPtsRef), rm.refPts.begin(),
         (dims_pts[0] * dims_pts[1] * sizeof(double)));
  memcpy((void *)PyArray_DATA(retPtsMov), rm.movPts.begin(),
         (dims_pts[0] * dims_pts[1] * sizeof(double)));

  /* "N" passes our references on to the tuple */
  return Py_BuildValue("NNN", retMat, retPtsRef, retPtsMov);
}

PyArrayObject *turbogreg_transform(PyObject *self, PyObject *args) {
  PyObject *mov;
  PyObject *mat;

  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "OO", &mov, &mat)) {
    return NULL;
  }

  /* Interpret the input objects as numpy arrays. */
  PyObject *po_mov_array = PyArray_FROM_OTF(mov, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (po_mov_array == NULL) {
    return NULL;
  }
  PyObject *po_mat_array = PyArray_FROM_OTF(mat, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (po_mov_array == NULL || po_mat_array == NULL) {
    Py_XDECREF(po_mat_array);
    Py_XDECREF(po_mov_array);
    return NULL;
  }

  if (!PyArray_Check(po_mat_array) || !PyArray_Check(po_mov_array)) {
    PyErr_SetString(PyExc_TypeError, "Expected a numpy array");
    return NULL;
  }

  PyArrayObject *mov_array = (PyArrayObject *)po_mov_array;
  PyArrayObject *mat_array = (PyArrayObject *)po_mat_array;

  int ndim_mov = (int)PyArray_NDIM(mov_array);
  int ndim_mat = (int)PyArray_NDIM(mat_array);

  if (ndim_mov != 2 || ndim_mat != 2) {
    Py_XDECREF(mat_array);
    Py_XDECREF(mov_array);
    PyErr_SetString(PyExc_ValueError, "Input arrays must be two dimensional");
    return NULL;
  }

  if (!checkImageSize(PyArray_DIM(mov_array, 0), PyArray_DIM(mov_array, 1), 1,
                      "The image to transform")) {
    Py_XDECREF(mat_array);
    Py_XDECREF(mov_array);
    return NULL;
  }

  /* How many data points are there? */
  int Nx_mov = (int)PyArray_DIM(mov_array, 0);
  int Ny_mov = (int)PyArray_DIM(mov_array, 1);
  int Nx_mat = (int)PyArray_DIM(mat_array, 0);
  int Ny_mat = (int)PyArray_DIM(mat_array, 1);

  if (Nx_mat != 2 || ((Ny_mat != 1) && (Ny_mat != 3) && (Ny_mat != 4))) {
    Py_XDECREF(mov_array);
    Py_XDECREF(mat_array);
    PyErr_SetString(
        PyExc_ValueError,
        "Transformation matrix must be of shape (2,1), (2,3) or (2,4)");
    return NULL;
  }

  /* Get pointers to the data as C-types. */
  double *img_mov = (double *)PyArray_DATA(mov_array);
  double *tmat = (double *)PyArray_DATA(mat_array);

  /* Reject matrices that are not finite or that map the image so far away
     that source coordinates would not fit into an int. */
  double bound[2];
  for (int i = 0; i < 2; i++) {
    const double *row = tmat + i * Ny_mat;
    bound[i] = std::abs(row[0]);
    if (Ny_mat >= 3) {
      bound[i] += std::abs(row[1]) * (double)Ny_mov +
                  std::abs(row[2]) * (double)Nx_mov;
    } else {
      bound[i] += (double)(i == 0 ? Ny_mov : Nx_mov);
    }
    if (Ny_mat == 4) {
      bound[i] += std::abs(row[3]) * (double)Ny_mov * (double)Nx_mov;
    }
  }
  /* the negated comparison also catches NaN */
  if (!(bound[0] < MAX_COORDINATE) || !(bound[1] < MAX_COORDINATE)) {
    Py_XDECREF(mov_array);
    Py_XDECREF(mat_array);
    PyErr_SetString(PyExc_ValueError,
                    "Transformation matrix must be finite and must not map "
                    "the image more than 1e9 pixels away");
    return NULL;
  }

  matrix<double> m(Nx_mat, Ny_mat);
  memcpy(m.begin(), tmat, (Nx_mat * Ny_mat * sizeof(double)));

  std::vector<double> imgout;
  NativeStatus status = NATIVE_OK;
  std::string message;

  Py_BEGIN_ALLOW_THREADS
  try {
    imgout =
        transformImg(m, img_mov, Ny_mov,
                     Nx_mov); // width and height (Nx/Ny) have to be swapped!
  } catch (const std::bad_alloc &) {
    status = NATIVE_NO_MEMORY;
  } catch (const std::exception &e) {
    status = NATIVE_ERROR;
    try {
      message = e.what();
    } catch (...) {
    }
  } catch (...) {
    status = NATIVE_ERROR;
  }
  Py_END_ALLOW_THREADS

  /* clean up */
  Py_XDECREF(mat_array);
  Py_XDECREF(mov_array);

  if (status != NATIVE_OK) {
    setNativeError(status, message);
    return NULL;
  }

  npy_intp dims[2] = {Nx_mov, Ny_mov};

  PyArrayObject *ret = (PyArrayObject*)PyArray_SimpleNew(2, (npy_intp *)&dims, NPY_DOUBLE);
  if (ret == NULL) {
    return NULL;
  }

  memcpy((void *)PyArray_DATA(ret), &imgout[0],
         (Nx_mov * Ny_mov * sizeof(double)));

  return ret;
}
