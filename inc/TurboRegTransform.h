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

#include "TurboReg.h"
#include "TurboRegImage.h"
#include "TurboRegMask.h"
#include "TurboRegPointHandler.h"
#include "matrix.h"
#include <string>



class TurboRegTransform
{

public:
    /*********************************************************************
     Keep a local copy of most everything. Select among the pre-stored
     constants.
     @param targetImg Target image pyramid.
     @param targetMsk Target mask pyramid.
     @param sourceImg Source image pyramid.
     @param sourceMsk Source mask pyramid.
     @param targetPh Target <code>TurboRegPointHandler</code> object.
     @param sourcePh Source <code>TurboRegPointHandler</code> object.
     @param transformation Transformation code.
     @param accelerated Trade-off between speed and accuracy.
     @param interactive Shows or hides the resulting image.
     ********************************************************************/
    TurboRegTransform (
		TurboRegImage *sourceImg,
		TurboRegMask *sourceMsk,
		TurboRegPointHandler *sourcePh,
		TurboRegImage *targetImg,
		TurboRegMask *targetMsk,
		TurboRegPointHandler *targetPh,
		Transformation transformation,
		bool accelerated
    );

    TurboRegTransform(
		TurboRegImage *sourceImg,
		TurboRegMask *sourceMsk,
		TurboRegPointHandler *sourcePh,
		Transformation transformation,
		bool accelerated
    );

    void appendTransformation (
        std::string pathAndFilename
    );

    void doBatchFinalTransform (
            std::vector<double> &pixels
    );
    std::vector<double> doFinalTransform (
            int width,
            int height
    );

	std::vector<double> doFinalTransform (
			TurboRegImage *sourceImg,
			TurboRegPointHandler *sourcePh,
			TurboRegImage *targetImg,
			TurboRegPointHandler *targetPh,
			Transformation transformation,
			bool accelerated
	);

	std::vector<double> doFinalTransform (
	        TurboRegImage *sourceImg,
			matrix<double> &m
	);

    matrix<double> getTransformationMatrix (
            matrix<double> &fromCoord,
            matrix<double> &toCoord
    );

    matrix<double> getTransformationMatrix ();

    void doRegistration (
    );

    void printMatrix(matrix<double> &m);
    void printPoints();


private:
/*....................................................................
	variables
....................................................................*/
    /*********************************************************************
     Maximal number of registration iterations per level, when
     speed is requested at the expense of accuracy. This number must be
     corrected so that there are more iterations at the coarse levels
     of the pyramid than at the fine levels.
     @see TurboRegTransform#ITERATION_PROGRESSION
     ********************************************************************/
	static const int FEW_ITERATIONS;

    /*********************************************************************
     Initial value of the Marquardt-Levenberg fudge factor.
     ********************************************************************/
	static const double FIRST_LAMBDA;

    /*********************************************************************
     Update parameter of the Marquardt-Levenberg fudge factor.
     ********************************************************************/
	static const double LAMBDA_MAGSTEP;

    /*********************************************************************
     Maximal number of registration iterations per level, when
     accuracy is requested at the expense of speed. This number must be
     corrected so that there are more iterations at the coarse levels
     of the pyramid than at the fine levels.
     @see TurboRegTransform#ITERATION_PROGRESSION
     ********************************************************************/
	static const int MANY_ITERATIONS;

    /*********************************************************************
     Minimal update distance of the landmarks, in pixel units, when
     accuracy is requested at the expense of speed. This distance does
     not depend on the pyramid level.
     ********************************************************************/
	static const double PIXEL_HIGH_PRECISION;

    /*********************************************************************
     Minimal update distance of the landmarks, in pixel units, when
     speed is requested at the expense of accuracy. This distance does
     not depend on the pyramid level.
     ********************************************************************/
	static const double PIXEL_LOW_PRECISION;

    /*********************************************************************
     Multiplicative factor that determines how many more iterations
     are allowed for a pyramid level one unit coarser.
     ********************************************************************/
	static const int ITERATION_PROGRESSION;

    bool accelerated;
    double c0;
    double c0u;
    double c0v;
    double c0uv;
    double c1;
    double c1u;
    double c1v;
    double c1uv;
    double c2;
    double c2u;
    double c2v;
    double c2uv;
    double c3;
    double c3u;
    double c3v;
    double c3uv;
    double pixelPrecision;
    double s;
    double t;
    double targetJacobian;
    double x;
    double y;
    matrix<double> sourcePoint;
    matrix<double> targetPoint;
    std::vector<double> dxWeight;
    std::vector<double> dyWeight;
    std::vector<double> xWeight;
    std::vector<double> yWeight;
    std::vector<int> xIndex;
    std::vector<int> yIndex;
    std::vector<double> inImg;
    std::vector<double> inMsk;
    std::vector<double> outImg;
    std::vector<double> outMsk;
    std::vector<double> xGradient;
    std::vector<double> yGradient;
    int inNx;
    int inNy;
    int iterationCost;
    int iterationPower;
    int maxIterations;
    int outNx;
    int outNy;
    int p;
    int pyramidDepth;
    int q;
    Transformation transformation;
    int twiceInNx;
    int twiceInNy;
    TurboRegImage *sourceImg;
    TurboRegImage *targetImg;
    TurboRegMask *sourceMsk;
    TurboRegMask *targetMsk;
    TurboRegPointHandler *sourcePh;
    bool bHasSourceMask;



    bool hasSourceMask() {
            return bHasSourceMask;
    }

    void affineTransform (
            matrix<double> &m
    );
    void affineTransform (
            matrix<double> &m,
            std::vector<double> &outMsk
    );
    void bilinearTransform (
            matrix<double> &m
    );
    void bilinearTransform (
            matrix<double> &m,
            std::vector<double> &outMsk
    );
    void computeBilinearGradientConstants (
    );
    double getAffineMeanSquares (
            matrix<double> &sourcePoint,
            matrix<double> &m
    );
    double getAffineMeanSquares (
            matrix<double> &sourcePoint,
            matrix<double> &m,
            std::vector<double> &gradient
    );
    double getAffineMeanSquares (
            matrix<double> &sourcePoint,
            matrix<double> &m,
            matrix<double> &hessian,
            std::vector<double> &gradient
    );
    double getBilinearMeanSquares (
            matrix<double> &m
    );
    double getBilinearMeanSquares (
            matrix<double> &m,
            matrix<double> &hessian,
            std::vector<double> &gradient
    );
    double getRigidBodyMeanSquares (
            matrix<double> &m
    );
    double getRigidBodyMeanSquares (
            matrix<double> &m,
            std::vector<double> &gradient

    );

    double getRigidBodyMeanSquares (
            matrix<double> &m,
            matrix<double> &hessian,
            std::vector<double> &gradient
    );
    double getScaledRotationMeanSquares (
            matrix<double> &sourcePoint,
            matrix<double> &m
    );

    double getScaledRotationMeanSquares (
            matrix<double> &sourcePoint,
            matrix<double> &m,
            std::vector<double> &gradient
    );

    double getScaledRotationMeanSquares (
            matrix<double> &sourcePoint,
            matrix<double> &m,
            matrix<double> &hessian,
            std::vector<double> &gradient
    );


    double getTranslationMeanSquares (
            matrix<double> &m
    );

    double getTranslationMeanSquares (
            matrix<double> &m,
            std::vector<double> &gradient
    );

    double getTranslationMeanSquares (
            matrix<double> &m,
            matrix<double> &hessian,
            std::vector<double> &gradient
    );

    double interpolate ();
    double interpolateDx ();
    double interpolateDy ();

    void inverseMarquardtLevenbergOptimization ();
    void inverseMarquardtLevenbergRigidBodyOptimization ();

    void invertGauss (
            matrix<double> &m
    );

    void marquardtLevenbergOptimization ();

    std::vector<double> matrixMultiply (
            matrix<double> &m,
            std::vector<double> &vector
    );

    void scaleBottomDownLandmarks ();
    void scaleUpLandmarks ();

    void translationTransform (
            matrix<double> &m
    );

    void translationTransform (
            matrix<double> &m,
            std::vector<double> &outMsk
    );

    void xDxWeights ();
    void xIndexes ();
    void xWeights ();
    void yDyWeights ();
    void yIndexes ();
    void yWeights ();

}; /* end class TurboRegTransform */
