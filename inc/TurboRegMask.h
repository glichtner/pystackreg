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

#include <stack>
#include "matrix.h"

#include "TurboRegImage.h"

#ifndef TURBOREGMASK_H_
#define TURBOREGMASK_H_

class MaskStackItem {
public:
    std::vector<double> halfMask;

    MaskStackItem(int size) {
        halfMask.resize(size);
    }
};

class TurboRegMask
{

public:
    TurboRegMask (TurboRegImage &img);
    TurboRegMask (matrix<double> &imp, int width, int height);

    void clearMask ();
    std::vector<double> &getMask();

    std::stack<MaskStackItem> getPyramid();

    void setPyramidDepth (int pyramidDepth);
    void init();

private:

    std::stack<MaskStackItem> pyramid;

    std::vector<double> mask;
    int width;
    int height;
    int pyramidDepth;

    void buildPyramid ();

    std::vector<double> getHalfMask2D (
            double *pFullMask,
            int fullWidth,
            int fullHeight,
            std::vector<double> &halfMask
    );


}; /* end class turboRegMask */

#endif
