/**
 * C++ Port of the TurboReg ImageJ Plugin
 * Original code by Philippe Thevenaz (see below)
 * Porting by Gregor Lichtner
 */


/*====================================================================
| Philippe Thevenaz
| EPFL/STI/IMT/LIB/BM.4.137
| Station 17
| CH-1015 Lausanne VD
| Switzerland
|
| phone (CET): +41(21)693.51.61
| fax: +41(21)693.37.01
| RFC-822: philippe.thevenaz@epfl.ch
| X-400: /C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
| URL: http://bigwww.epfl.ch/
\===================================================================*/

/*====================================================================
| This work is based on the following paper:
|
| P. Thevenaz, U.E. Ruttimann, M. Unser
| A Pyramid Approach to Subpixel Registration Based on Intensity
| IEEE Transactions on Image Processing
| vol. 7, no. 1, pp. 27-41, January 1998.
|
| This paper is available on-line at
| http://bigwww.epfl.ch/publications/thevenaz9801.html
|
| Other relevant on-line publications are available at
| http://bigwww.epfl.ch/publications/
\===================================================================*/

/*====================================================================
| Additional help available at http://bigwww.epfl.ch/thevenaz/turboreg/
|
| You'll be free to use this software for research purposes, but you
| should not redistribute it without our consent. In addition, we expect
| you to include a citation or acknowledgment whenever you present or
| publish results that are based on it.
\===================================================================*/

#include <stdexcept>

#include "TurboReg.h"

int getPyramidDepth (
        int sw,
        int sh,
        int tw,
        int th
) {
    int pyramidDepth = 1;
    while (((2 * PYRAMID_MIN_SIZE) <= sw)
            && ((2 * PYRAMID_MIN_SIZE) <= sh)
            && ((2 * PYRAMID_MIN_SIZE) <= tw)
            && ((2 * PYRAMID_MIN_SIZE) <= th)) {
        sw /= 2;
        sh /= 2;
        tw /= 2;
        th /= 2;
        pyramidDepth++;
    }
    return(pyramidDepth);
} /* end getPyramidDepth */

Transformation getTransformationFromMatrix(matrix<double> &m) {

	Transformation transformation;

	switch(m.ncols()) {
    case 1: transformation = TRANSLATION; break;
    case 3: transformation = AFFINE; break; // or RIGID_BODY or SCALED_ROT, but that doesn't matter
    case 4: transformation = BILINEAR; break;
    default:
    	throw std::runtime_error("Invalid transformation");

    }

	return transformation;
}
