#include <cmath>
#include "TurboRegPointHandler.h"

#define TURBOREG_MODE

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
	precisionPoint{precisionPoint}
{

};


void TurboRegPointHandler::setPointsByTransformation(
		int width,
		int height,
		Transformation transformation
) {

	switch (transformation) {
	case TRANSLATION: {	//AFFINE: { //should be three points, not one.
		precisionPoint.resize(1,2);
		precisionPoint(0, 0) = (double)(width / 2);
		precisionPoint(0, 1) = (double)(height / 2);
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
		precisionPoint(0, 0) = (double)(width / 2);
		precisionPoint(0, 1) = (double)(height / 2);
		precisionPoint(1, 0) = (double)(width / 2);
		precisionPoint(1, 1) = (double)(height / 4);
		precisionPoint(2, 0) = (double)(width / 2);
		precisionPoint(2, 1) = (double)((3 * height) / 4);
#endif
		break;

	}
	case SCALED_ROTATION: { //two points
		precisionPoint.resize(2,2);
		precisionPoint(0, 0) = (double)(width / 4);
		precisionPoint(0, 1) = (double)(height / 2);
		precisionPoint(1, 0) = (double)((3 * width) / 4);
		precisionPoint(1, 1) = (double)(height / 2);
		break;
	}
	case AFFINE: {
		precisionPoint.resize(3,2);
		precisionPoint(0, 0) = (double)(width / 2);
		precisionPoint(0, 1) = (double)(height / 4);
		precisionPoint(1, 0) = (double)(width / 4);
		precisionPoint(1, 1) = (double)((3 * height) / 4);
		precisionPoint(2, 0) = (double)((3 * width) / 4);
		precisionPoint(2, 1) = (double)((3 * height) / 4);
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
