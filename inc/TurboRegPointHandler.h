#include "TurboReg.h"
#include "TurboRegImage.h"
#include <cmath>
#include "matrix.h"

#ifndef TURBOREGPOINTHANDLER_H_
#define TURBOREGPOINTHANDLER_H_

class TurboRegPointHandler
{ /* class TurboRegPointHandler */
public:

    TurboRegPointHandler (
        TurboRegImage &img,
		Transformation transformation
    );

    TurboRegPointHandler (
    	matrix<double> &precisionPoint
    );

    matrix<double> &getPoints();

    void setPoints (
        matrix<double> &precisionPoint
    );

    void setPointsByTransformation(
		int width,
		int height,
    	Transformation transformation
	);

private:
    /*********************************************************************
     The number of points we are willing to deal with is at most
    <code>4</code>.
    @see TurboRegDialog#transformation
    ********************************************************************/
    int NUM_POINTS = 4;

    /*********************************************************************
     The golden ratio mathematical constant determines where to put the
    initial landmarks.
    ********************************************************************/
    double GOLDEN_RATIO = 0.5 * (sqrt(5.0) - 1.0);


    matrix<double> precisionPoint; // = new double[NUM_POINTS * 2];
};


#endif
