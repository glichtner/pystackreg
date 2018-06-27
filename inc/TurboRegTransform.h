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
            TurboRegImage &sourceImg,
            TurboRegMask &sourceMsk,
            TurboRegPointHandler &sourcePh,
            TurboRegImage &targetImg,
            TurboRegMask &targetMsk,
            TurboRegPointHandler &targetPh,
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
			TurboRegImage &sourceImg,
			TurboRegPointHandler &sourcePh,
			TurboRegImage &targetImg,
			TurboRegPointHandler &targetPh,
			Transformation transformation,
			bool accelerated
	);

	std::vector<double> doFinalTransform (
	        TurboRegImage &sourceImg,
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
    const int FEW_ITERATIONS = 5;

    /*********************************************************************
     Initial value of the Marquardt-Levenberg fudge factor.
     ********************************************************************/
    const double FIRST_LAMBDA = 1.0;

    /*********************************************************************
     Update parameter of the Marquardt-Levenberg fudge factor.
     ********************************************************************/
    const double LAMBDA_MAGSTEP = 4.0;

    /*********************************************************************
     Maximal number of registration iterations per level, when
     accuracy is requested at the expense of speed. This number must be
     corrected so that there are more iterations at the coarse levels
     of the pyramid than at the fine levels.
     @see TurboRegTransform#ITERATION_PROGRESSION
     ********************************************************************/
    const int MANY_ITERATIONS = 10;

    /*********************************************************************
     Minimal update distance of the landmarks, in pixel units, when
     accuracy is requested at the expense of speed. This distance does
     not depend on the pyramid level.
     ********************************************************************/
    const double PIXEL_HIGH_PRECISION = 0.001;

    /*********************************************************************
     Minimal update distance of the landmarks, in pixel units, when
     speed is requested at the expense of accuracy. This distance does
     not depend on the pyramid level.
     ********************************************************************/
    const double PIXEL_LOW_PRECISION = 0.1;

    /*********************************************************************
     Multiplicative factor that determines how many more iterations
     are allowed for a pyramid level one unit coarser.
     ********************************************************************/
    const int ITERATION_PROGRESSION = 2;

    bool accelerated;
    double c0 = 0;
    double c0u = 0;
    double c0v = 0;
    double c0uv = 0;
    double c1 = 0;
    double c1u = 0;
    double c1v = 0;
    double c1uv = 0;
    double c2 = 0;
    double c2u = 0;
    double c2v = 0;
    double c2uv = 0;
    double c3 = 0;
    double c3u = 0;
    double c3v = 0;
    double c3uv = 0;
    double pixelPrecision = 0;
    double s = 0;
    double t = 0;
    double targetJacobian = 0;
    double x = 0;
    double y = 0;
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
    TurboRegImage &sourceImg;
    TurboRegImage &targetImg;
    TurboRegMask &sourceMsk;
    TurboRegMask &targetMsk;
    TurboRegPointHandler &sourcePh;
    bool bHasSourceMask = true;



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
