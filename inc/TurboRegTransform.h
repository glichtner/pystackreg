#include "TurboRegImage.h"
#include "TurboRegMask.h"
#include "TurboRegPointHandler.h"


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
            TurboRegImage sourceImg,
            TurboRegMask sourceMsk,
            TurboRegPointHandler sourcePh,
            TurboRegImage targetImg,
            TurboRegMask targetMsk,
            TurboRegPointHandler targetPh,
            int transformation,
            bool accelerated,
            bool interactive
    );

    appendTransformation (
        final String pathAndFilename
    );

    void doBatchFinalTransform (
            double[] pixels
    );
    ImagePlus doFinalTransform (
            int width,
            int height
    );
    double[] doFinalTransform (
            TurboRegImage sourceImg,
            TurboRegPointHandler sourcePh,
            TurboRegImage targetImg,
            TurboRegPointHandler targetPh,
            int transformation,
            bool accelerated
    );
    void doRegistration (
    );

 

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
    bool interactive;
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
    double[][] sourcePoint;
    double[][] targetPoint;
    double* dxWeight = new double[4];
    double* dyWeight = new double[4];
    double* xWeight = new double[4];
    double* yWeight = new double[4];
    int* xIndex = new int[4];
    int* yIndex = new int[4];
    double* inImg;
    double* inMsk;
    double* outImg;
    double* outMsk;
    double* xGradient;
    double* yGradient;
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
    int transformation;
    int twiceInNx;
    int twiceInNy;
    TurboRegImage sourceImg;
    TurboRegImage targetImg;
    TurboRegMask sourceMsk;
    TurboRegMask targetMsk;
    TurboRegPointHandler sourcePh;


    void affineTransform (
            double[][] matrix
    );
    void affineTransform (
            double[][] matrix,
            double[] outMsk
    );
    void bilinearTransform (
            double[][] matrix
    );
    void bilinearTransform (
            double[][] matrix,
            double[] outMsk
    );
    void computeBilinearGradientConstants (
    );
    double getAffineMeanSquares (
            double[][] sourcePoint,
            double[][] matrix
    );
    double getAffineMeanSquares (
            double[][] sourcePoint,
            double[][] matrix,
            double[] gradient
    );
    double getAffineMeanSquares (
            double[][] sourcePoint,
            double[][] matrix,
            double[][] hessian,
            double[] gradient
    );
    double getBilinearMeanSquares (
            double[][] matrix
    );
    double getBilinearMeanSquares (
            double[][] matrix,
            double[][] hessian,
            double[] gradient
    );
    double getRigidBodyMeanSquares (
            double[][] matrix
    );
    double getRigidBodyMeanSquares (
            double[][] matrix,
            double[] gradient

    );

    double getRigidBodyMeanSquares (
            double[][] matrix,
            double[][] hessian,
            double[] gradient
    );
    double getScaledRotationMeanSquares (
            double[][] sourcePoint,
            double[][] matrix
    );

    double getScaledRotationMeanSquares (
            double[][] sourcePoint,
            double[][] matrix,
            double[] gradient
    );

    double getScaledRotationMeanSquares (
            double[][] sourcePoint,
            double[][] matrix,
            double[][] hessian,
            double[] gradient
    );

    double[][] getTransformationMatrix (
            double[][] fromCoord,
            double[][] toCoord
    );

    double getTranslationMeanSquares (
            double[][] matrix
    );

    double getTranslationMeanSquares (
            double[][] matrix,
            double[] gradient
    );

    double getTranslationMeanSquares (
            double[][] matrix,
            double[][] hessian,
            double[] gradient
    );

    double interpolate (
    );

    double interpolateDx (
    );

    double interpolateDy (
    );

    void inverseMarquardtLevenbergOptimization (
            int workload
    );

    void inverseMarquardtLevenbergRigidBodyOptimization (
            int workload
    );

    void invertGauss (
            double[][] matrix
    );

    void marquardtLevenbergOptimization (
            int workload
    );
    double[] matrixMultiply (
            double[][] matrix,
            double[] vector
    );

    void scaleBottomDownLandmarks (
    );

    void scaleUpLandmarks (
    );

    void translationTransform (
            double[][] matrix
    );

    void translationTransform (
            double[][] matrix,
            double[] outMsk
    );

    void xDxWeights (
    );

    void xIndexes (
    );

    void xWeights (
    );

    void yDyWeights (
    );

    void yIndexes (
    );

    void yWeights (
    );

} /* end class TurboRegTransform */