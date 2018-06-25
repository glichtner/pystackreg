#include "TurboRegImage.h"
#include <cmath>

class TurboRegPointHandler
{ /* class TurboRegPointHandler */

public:
    TurboRegPointHandler (
        TurboRegImage img
    );

    double* getPoints ();

    void setPoints (
        double* precisionPoint
    );


private:

    /*....................................................................
        variables
    ....................................................................*/
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



    double* precisionPoint = new double[NUM_POINTS * 2];
    int transformation;
    int currentPoint = 0;

}
