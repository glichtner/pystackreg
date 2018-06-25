#include "TurboRegPointHandler.h"




TurboRegPointHandler::TurboRegPointHandler (
	TurboRegImage img
) {
    int width = img.getWidth();
    int height = img.getHeight();
}


double* TurboRegPointHandler::getPoints (
) {
	if (interactive) {
		if (transformation == turboRegDialog.RIGID_BODY) {
			double[][] points = new double[transformation][2];
			for (int k = 0; (k < transformation); k++) {
				points[k][0] = (double)point[k].x;
				points[k][1] = (double)point[k].y;
			}
			return(points);
		}
		else {
			double[][] points = new double[transformation / 2][2];
			for (int k = 0; (k < (transformation / 2)); k++) {
				points[k][0] = (double)point[k].x;
				points[k][1] = (double)point[k].y;
			}
			return(points);
		}
	}
	else {
		return(precisionPoint);
	}
} /* end getPoints */

void TurboRegPointHandler::setPoints (
	double[][] precisionPoint
) {
	interactive = false;
	if (transformation == turboRegDialog.RIGID_BODY) {
		for (int k = 0; (k < transformation); k++) {
			point[k].x = (int)Math.round(precisionPoint[k][0]);
			point[k].y = (int)Math.round(precisionPoint[k][1]);
			this.precisionPoint[k][0] = precisionPoint[k][0];
			this.precisionPoint[k][1] = precisionPoint[k][1];
		}
	}
	else {
		for (int k = 0; (k < (transformation / 2)); k++) {
			point[k].x = (int)Math.round(precisionPoint[k][0]);
			point[k].y = (int)Math.round(precisionPoint[k][1]);
			this.precisionPoint[k][0] = precisionPoint[k][0];
			this.precisionPoint[k][1] = precisionPoint[k][1];
		}
	}
} /* end setPoints */
