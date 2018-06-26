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
        TurboRegImage img
    );

    matrix<double> &getPoints();

    void setPoints (
        matrix<double> &precisionPoint
    );


private:
    /*********************************************************************
     The magnifying tool is set in eleventh position to be coherent with
    ImageJ.
    ********************************************************************/
    const int MAGNIFIER = 11;

    /*********************************************************************
     The moving tool is set in second position to be coherent with the
    <code>PointPicker_</code> plugin.
    ********************************************************************/
    const int MOVE_CROSS = 1;

    /*********************************************************************
     The number of points we are willing to deal with is at most
    <code>4</code>.
    @see TurboRegDialog#transformation
    ********************************************************************/
    const int NUM_POINTS = 4;

    /*....................................................................
        variables
    ....................................................................*/
    /*********************************************************************
     Serialization version number.
    ********************************************************************/
    const long serialVersionUID = 1L;

    /*********************************************************************
     The drawn landmarks fit in a 11x11 matrix.
    ********************************************************************/
    const int CROSS_HALFSIZE = 5;

    /*********************************************************************
     The golden ratio mathematical constant determines where to put the
    initial landmarks.
    ********************************************************************/
    const double GOLDEN_RATIO = 0.5 * (sqrt(5.0) - 1.0);

    matrix<double> precisionPoint; // = new double[NUM_POINTS * 2];
};


#endif
