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


#ifndef TURBOREG_H_
#define TURBOREG_H_

#include "matrix.h"

typedef enum Transformation {
/*********************************************************************
 Generic geometric transformation.
 ********************************************************************/
//GENERIC_TRANSFORMATION = -1,

/*********************************************************************
 A translation is described by a single point. It keeps area, angle,
 and orientation. A translation is determined by 2 parameters.
 ********************************************************************/
TRANSLATION = 2,

/*********************************************************************
 A single points determines the translation component of a rigid-body
 transformation. As the rotation is given by a scalar number, it is
 not natural to represent this number by coordinates of a point. The
 rigid-body transformation is determined by 3 parameters.
 ********************************************************************/
RIGID_BODY = 3,

/*********************************************************************
 A pair of points determines the combination of a translation, of
 a rotation, and of an isotropic scaling. Angles are conserved. A
 scaled rotation is determined by 4 parameters.
 ********************************************************************/
SCALED_ROTATION = 4,

/*********************************************************************
 Three points generate an affine transformation, which is any
 combination of translation, rotation, isotropic scaling, anisotropic
 scaling, shearing, and skewing. An affine transformation maps
 parallel lines onto parallel lines and is determined by 6 parameters.
 ********************************************************************/
AFFINE = 6,


/*********************************************************************
 Four points describe a bilinear transformation, where a point of
 coordinates (x, y) is mapped on a point of coordinates (u, v) such
 that u = p0 + p1 x + p2 y + p3 x y and v = q0 + q1 x + q2 y + q3 x y.
 Thus, u and v are both linear in x, and in y as well. The bilinear
 transformation is determined by 8 parameters.
 ********************************************************************/
BILINEAR = 8,
} Transformation;

/*********************************************************************
 Minimal linear dimension of an image in the multiresolution pyramid.
 ********************************************************************/
const static int PYRAMID_MIN_SIZE = 12;


int getPyramidDepth(int sourceWidth, int sourceHeight, int targetWidth, int targetHeight);
Transformation getTransformationFromMatrix(matrix<double> &m);

#endif
