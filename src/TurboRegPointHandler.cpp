#include <cmath>
#include "TurboRegPointHandler.h"

TurboRegPointHandler::TurboRegPointHandler (
	TurboRegImage img
) : precisionPoint(NUM_POINTS, 2)

{
    int width = img.getWidth();
    int height = img.getHeight();

	precisionPoint(0,0) = floor(width / 2.0);
	precisionPoint(0,1) = floor(height / 2.0);
	
	precisionPoint(1,0) = floor(width / 2.0);
	precisionPoint(1,1) = ceil( height * 0.25 * GOLDEN_RATIO);

	precisionPoint(2,0) = floor(width / 2.0);
	precisionPoint(2,1) = height - ceil(height * 0.25 * GOLDEN_RATIO);

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
