/*
 * PyStackReg.cpp
 *
 *  Created on: Jun 27, 2018
 *      Author: lichtneg
 */


#include <Python.h>
#include <numpy/arrayobject.h>


#include <iostream>
#include <vector>
#include <cstdint>
#include <iterator>

#include "TurboReg.h"
#include "TurboRegImage.h"
#include "TurboRegMask.h"
#include "TurboRegPointHandler.h"
#include "TurboRegTransform.h"
#include "matrix.h"


class PyStackReg {

public:
	int val;
	PyStackReg() :
		transformation(RIGID_BODY),
		val(42) {}


	PyStackReg(Transformation transformation)
: transformation{transformation}
	{}

	void registerImage(double *pDataRef, double *pDataMov, int width, int height) {
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

		isregistered = true;
		tmat = tform.getTransformationMatrix();

		/*std::vector<double> imgout = tform.doFinalTransform(width, height);

		TurboRegPointHandler refPH2(refPH.getPoints());
		TurboRegPointHandler movPH2(movPH.getPoints());

		std::vector<double> imgout2 = tform.doFinalTransform (movImg, movPH2, refImg, refPH2, RIGID_BODY, false);

		matrix<double> tm = tform.getTransformationMatrix();*/

	}

private:
	Transformation transformation;
	bool isregistered = false;
	matrix<double> tmat;

};





#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    void destroy_funkster(PyObject* funkster) {
        delete (PyStackReg*)PyCapsule_GetPointer(funkster, "val");
    }

    static PyObject* get_funkster_val(PyObject *self, PyObject *args) {
        PyObject* pf = NULL;
        if(!PyArg_UnpackTuple(args, "funk", 1, 1, &pf)) return NULL;
        PyStackReg* f = (PyStackReg*)PyCapsule_GetPointer(pf, "val");
        return PyLong_FromLong(f->val);
    }

    static PyObject* make_funkster(PyObject *self, PyObject *args) {
        return PyCapsule_New((void*)new PyStackReg(), "val", destroy_funkster);
    }

    static PyMethodDef module_methods[] = {
        {"make_funkster", make_funkster, METH_NOARGS, "Make a Funkster."},
        {"get_funkster_val", get_funkster_val, METH_VARARGS, "Get val from Funkster"},
        {NULL, NULL, 0, NULL} /* Sentinel */
    };

    static struct PyModuleDef pystackreg =
    {
        PyModuleDef_HEAD_INIT,
        "pystackreg", /* name of module */
        "pystackreg\n", /* module documentation, may be NULL */
        -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        module_methods
    };

    PyMODINIT_FUNC PyInit_pystackreg(void)
    {
        /* Load `numpy` functionality. */
        import_array();

        return PyModule_Create(&pystackreg);
    }
#ifdef __cplusplus
} // end extern "C"
#endif // __cplusplus



