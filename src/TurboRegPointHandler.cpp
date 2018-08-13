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

#include <cmath>
#include "TurboRegPointHandler.h"

//#define TURBOREG_MODE

double TurboRegPointHandler::GOLDEN_RATIO = 0.5 * (sqrt(5.0) - 1.0);

TurboRegPointHandler::TurboRegPointHandler (
	TurboRegImage &img,
	Transformation transformation
) :
	precisionPoint()
{
	setPointsByTransformation(img.getWidth(), img.getHeight(), transformation);
};

TurboRegPointHandler::TurboRegPointHandler (
	matrix<double> &precisionPoint
) :
	precisionPoint(precisionPoint)
{

};


void TurboRegPointHandler::setPointsByTransformation(
		int width,
		int height,
		Transformation transformation
) {
	// increment width and height because images are supplied cropped by 1 px from python
	// (this cropping is exactly what StackReg does when calling TurboReg)
	width++;
	height++;

	switch (transformation) {
	case TRANSLATION: {	//AFFINE: { //should be three points, not one.
		precisionPoint.resize(1,2);
		precisionPoint(0, 0) = floor((double)(width / 2.0));
		precisionPoint(0, 1) = floor((double)(height / 2.0));
		break;
	}
	case RIGID_BODY: { //three points

		precisionPoint.resize(3,2);

#ifdef TURBOREG_MODE
		precisionPoint(0,0) = floor(width / 2.0);
		precisionPoint(0,1) = floor(height / 2.0);
		precisionPoint(1,0) = floor(width / 2.0);
		precisionPoint(1,1) = ceil( height * 0.25 * GOLDEN_RATIO);
		precisionPoint(2,0) = floor(width / 2.0);
		precisionPoint(2,1) = height - ceil(height * 0.25 * GOLDEN_RATIO);
#else
		precisionPoint(0, 0) = floor((double)(width / 2.0));
		precisionPoint(0, 1) = floor((double)(height / 2.0));
		precisionPoint(1, 0) = floor((double)(width / 2.0));
		precisionPoint(1, 1) = floor((double)(height / 4.0));
		precisionPoint(2, 0) = floor((double)(width / 2.0));
		precisionPoint(2, 1) = floor((double)((3.0 * height) / 4.0));
#endif
		break;

	}
	case SCALED_ROTATION: { //two points
		precisionPoint.resize(2,2);
		precisionPoint(0, 0) = floor((double)(width / 4.0));
		precisionPoint(0, 1) = floor((double)(height / 2.0));
		precisionPoint(1, 0) = floor((double)((3.0 * width) / 4.0));
		precisionPoint(1, 1) = floor((double)(height / 2.0));
		break;
	}
	case AFFINE: {
		precisionPoint.resize(3,2);
#ifdef TURBOREG_MODE
		precisionPoint(0,0) = floor(width / 2.0);
		precisionPoint(0,1) = floor(0.25 * GOLDEN_RATIO * (double)height);
		precisionPoint(1,0) = floor(0.25 * GOLDEN_RATIO * (double)width);
		precisionPoint(1,1) = height - ceil(0.25 * GOLDEN_RATIO * height);
		precisionPoint(2,0) = width - ceil(0.25 * GOLDEN_RATIO * (double)width);
		precisionPoint(2,1) = height - ceil(height * 0.25 * GOLDEN_RATIO);

#else
		precisionPoint(0, 0) = floor((double)(width / 2.0));
		precisionPoint(0, 1) = floor((double)(height / 4.0));
		precisionPoint(1, 0) = floor((double)(width / 4.0));
		precisionPoint(1, 1) = floor((double)((3.0 * height) / 4.0));
		precisionPoint(2, 0) = floor((double)((3.0 * width) / 4.0));
		precisionPoint(2, 1) = floor((double)((3.0 * height) / 4.0));

#endif
		break;
	}

	case BILINEAR: {
		precisionPoint.resize(4,2);

		precisionPoint(0, 0) = (floor(0.25 * GOLDEN_RATIO * (double)width));
		precisionPoint(0, 1) = (floor(0.25 * GOLDEN_RATIO * (double)height));
		precisionPoint(1, 0) = (floor(0.25 * GOLDEN_RATIO * (double)width));
		precisionPoint(1, 1) = height - (ceil(0.25 * GOLDEN_RATIO * (double)height));
		precisionPoint(2, 0) = width - (ceil(0.25 * GOLDEN_RATIO * (double)width));
		precisionPoint(2, 1) = (floor(0.25 * GOLDEN_RATIO * (double)height));
		precisionPoint(3, 0) = width - (ceil(0.25 * GOLDEN_RATIO * (double)width));
		precisionPoint(3, 1) = height - (ceil(0.25 * GOLDEN_RATIO * (double)height));

		break;
	}

	default: {

		return;
	}

	}
}

matrix<double> &TurboRegPointHandler::getPoints (
) {
	return(precisionPoint);
} /* end getPoints */

void TurboRegPointHandler::setPoints (
	matrix<double> &precisionPoint
) {

	this->precisionPoint = precisionPoint;
	/*if (transformation == RIGID_BODY) {
		for (int k = 0; (k < transformation); k++) {
			point[k].x = (int)Math.round(precisionPoint[k][0]);
			point[k].y = (int)Math.round(precisionPoint[k][1]);
			this->precisionPoint[k][0] = precisionPoint[k][0];
			this->precisionPoint[k][1] = precisionPoint[k][1];
		}
	}
	else {
		for (int k = 0; (k < (transformation / 2)); k++) {
			point[k].x = (int)Math.round(precisionPoint[k][0]);
			point[k].y = (int)Math.round(precisionPoint[k][1]);
			this->precisionPoint[k][0] = precisionPoint[k][0];
			this->precisionPoint[k][1] = precisionPoint[k][1];
		}
	}*/
} /* end setPoints */
