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
#define _SECURE_SCL_DEPRECATE 0
#pragma warning(disable:4996)
#include <cmath>
#include <algorithm>
#include "TurboRegTransform.h"
#include "matrix.h"

#ifdef _MSC_VER

template <typename T, size_t N>
T* begin(T(&arr)[N]) { return &arr[0]; }
template <typename T, size_t N>
T* end(T(&arr)[N]) { return &arr[0] + N; }
#endif

const int TurboRegTransform::FEW_ITERATIONS = 5;
const double TurboRegTransform::FIRST_LAMBDA = 1.0;
const double TurboRegTransform::LAMBDA_MAGSTEP = 4.0;
const int TurboRegTransform::MANY_ITERATIONS = 10;
const double TurboRegTransform::PIXEL_HIGH_PRECISION = 0.001;
const double TurboRegTransform::PIXEL_LOW_PRECISION = 0.1;
const int TurboRegTransform::ITERATION_PROGRESSION = 2;

/*....................................................................
constructors
....................................................................*/
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
TurboRegTransform::TurboRegTransform(
        TurboRegImage *sourceImg,
        TurboRegMask *sourceMsk,
        TurboRegPointHandler *sourcePh,
        TurboRegImage *targetImg,
        TurboRegMask *targetMsk,
        TurboRegPointHandler *targetPh,
        Transformation transformation,
        bool accelerated
) :
	accelerated(accelerated),
	dxWeight(4),
    dyWeight(4),
    xWeight(4),
    yWeight(4),
    xIndex(4),
    yIndex(4),
	transformation(transformation),
    sourceImg(sourceImg),
	targetImg(targetImg),
	sourceMsk(sourceMsk),
	targetMsk(targetMsk),
	sourcePh(sourcePh),
	bHasSourceMask(true)
 {
    sourcePoint = sourcePh->getPoints();
    targetPoint = targetPh->getPoints();

    if (accelerated) {
        pixelPrecision = PIXEL_LOW_PRECISION;
        maxIterations = FEW_ITERATIONS;
    }
    else {
        pixelPrecision = PIXEL_HIGH_PRECISION;
        maxIterations = MANY_ITERATIONS;
    }
} /* end TurboRegTransform */


 TurboRegTransform::TurboRegTransform(
         TurboRegImage *sourceImg,
         TurboRegMask *sourceMsk,
         TurboRegPointHandler *sourcePh,
		 Transformation transformation,
		 bool accelerated
 ) :
	 accelerated(accelerated),
     dxWeight(4),
     dyWeight(4),
     xWeight(4),
     yWeight(4),
     xIndex(4),
     yIndex(4),
	 transformation(transformation),
	 sourceImg(sourceImg),
	 sourceMsk(sourceMsk),
	 sourcePh(sourcePh),
	 bHasSourceMask(true)
  {
     pixelPrecision = PIXEL_HIGH_PRECISION;
     maxIterations = MANY_ITERATIONS;

 } /* end TurboRegTransform */

/*....................................................................
public methods
....................................................................*/
/*********************************************************************
 Append the current landmarks into a text file. Rigid format.
    @param pathAndFilename Path and name of the file where batch results
    are being written.
    @see TurboRegDialog#loadLandmarks()
    ********************************************************************/
void TurboRegTransform::appendTransformation (
        std::string pathAndFilename
) {
    outNx = targetImg->getWidth();
    outNy = targetImg->getHeight();
    inNx = sourceImg->getWidth();
    inNy = sourceImg->getHeight();
    /*if (pathAndFilename == null) {
        return;
    }
    FileWriter fw = new FileWriter(pathAndFilename, true);
    fw.write("\n");
    switch (transformation) {
        case TRANSLATION: {
            fw.write("TRANSLATION\n");
            break;
        }
        case RIGID_BODY: {
            fw.write("RIGID_BODY\n");
            break;
        }
        case SCALED_ROTATION: {
            fw.write("SCALED_ROTATION\n");
            break;
        }
        case AFFINE: {
            fw.write("AFFINE\n");
            break;
        }
        case BILINEAR: {
            fw.write("BILINEAR\n");
            break;
        }
    }
    fw.write("\n");
    fw.write("Source size\n");
    fw.write(inNx + "\t" + inNy + "\n");
    fw.write("\n");
    fw.write("Target size\n");
    fw.write(outNx + "\t" + outNy + "\n");
    fw.write("\n");
    fw.write("Refined source landmarks\n");
    if (transformation == RIGID_BODY) {
        for (int i = 0; (i < transformation); i++) {
            fw.write(sourcePoint(i, 0) + "\t" + sourcePoint(i, 1) + "\n");
        }
    }
    else {
        for (int i = 0; (i < (transformation / 2)); i++) {
            fw.write(sourcePoint(i, 0) + "\t" + sourcePoint(i, 1) + "\n");
        }
    }
    fw.write("\n");
    fw.write("Target landmarks\n");
    if (transformation == RIGID_BODY) {
        for (int i = 0; (i < transformation); i++) {
            fw.write(targetPoint(i, 0) + "\t" + targetPoint(i, 1) + "\n");
        }
    }
    else {
        for (int i = 0; (i < (transformation / 2)); i++) {
            fw.write(targetPoint(i, 0) + "\t" + targetPoint(i, 1) + "\n");
        }
    }
    fw.close();*/

} /* end appendTransformation */

/*********************************************************************
 Compute the image.
    ********************************************************************/
void TurboRegTransform::doBatchFinalTransform (
        std::vector<double> &pixels
) {
    if (accelerated) {
        inImg = sourceImg->getImage();
    }
    else {
        inImg = sourceImg->getCoefficient();
    }
    inNx = sourceImg->getWidth();
    inNy = sourceImg->getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;
    outImg = pixels;
    outNx = targetImg->getWidth();
    outNy = targetImg->getHeight();
    matrix<double> m = getTransformationMatrix(targetPoint, sourcePoint);
    switch (transformation) {
        case TRANSLATION: {
            translationTransform(m);
            break;
        }
        case RIGID_BODY:
        case SCALED_ROTATION:
        case AFFINE: {
            affineTransform(m);
            break;
        }
        case BILINEAR: {
            bilinearTransform(m);
            break;
        }
    }
} /* end doBatchFinalTransform */

/*********************************************************************
 Compute the image.
    ********************************************************************/
std::vector<double> TurboRegTransform::doFinalTransform (
        int width,
        int height
) {
    if (accelerated) {
        inImg = sourceImg->getImage();
    }
    else {
        inImg = sourceImg->getCoefficient();
    }
    inMsk = sourceMsk->getMask();
    inNx = sourceImg->getWidth();
    inNy = sourceImg->getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;

    outImg.resize(width * height, 0);
    outMsk.resize(width * height, 1);

    outNx = width;
    outNy = height;

    matrix<double> m = getTransformationMatrix(targetPoint, sourcePoint);

    switch (transformation) {
        case TRANSLATION: {
            translationTransform(m, outMsk);
            break;
        }
        case RIGID_BODY:
        case SCALED_ROTATION:
        case AFFINE: {
            affineTransform(m, outMsk);
            break;
        }
        case BILINEAR: {
            bilinearTransform(m, outMsk);
            break;
        }
    }

    return(outImg);
} /* end doFinalTransform */

/*********************************************************************
 Compute the image.
    ********************************************************************/
std::vector<double> TurboRegTransform::doFinalTransform (
        TurboRegImage *sourceImg,
        TurboRegPointHandler *sourcePh,
        TurboRegImage *targetImg,
        TurboRegPointHandler *targetPh,
        Transformation transformation,
        bool accelerated
) {
    this->sourceImg = sourceImg;
    this->targetImg = targetImg;
    this->sourcePh = sourcePh;
    this->transformation = transformation;
    this->accelerated = accelerated;
    sourcePoint = sourcePh->getPoints();
    targetPoint = targetPh->getPoints();
    if (accelerated) {
        inImg = sourceImg->getImage();
    }
    else {
        inImg = sourceImg->getCoefficient();
    }
    inNx = sourceImg->getWidth();
    inNy = sourceImg->getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;
    outNx = targetImg->getWidth();
    outNy = targetImg->getHeight();
    outImg.resize(outNx * outNy);

    matrix<double> m = getTransformationMatrix(targetPoint, sourcePoint);

    switch (transformation) {
        case TRANSLATION: {
            translationTransform(m);
            break;
        }
        case RIGID_BODY:
        case SCALED_ROTATION:
        case AFFINE: {
            affineTransform(m);
            break;
        }
        case BILINEAR: {
            bilinearTransform(m);
            break;
        }
    }
    return(outImg);
} /* end doFinalTransform */


std::vector<double> TurboRegTransform::doFinalTransform (
        TurboRegImage *sourceImg,
		matrix<double> &m
) {
    this->sourceImg = sourceImg;
    if (accelerated) {
        inImg = sourceImg->getImage();
    }
    else {
        inImg = sourceImg->getCoefficient();
    }
    inNx = sourceImg->getWidth();
    inNy = sourceImg->getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;

    outNx = sourceImg->getWidth();
    outNy = sourceImg->getHeight();
    outImg.resize(outNx * outNy);

    switch (transformation) {
        case TRANSLATION: {
            translationTransform(m);
            break;
        }
        case RIGID_BODY:
        case SCALED_ROTATION:
        case AFFINE: {
            affineTransform(m);
            break;
        }
        case BILINEAR: {
            bilinearTransform(m);
            break;
        }
    }
    return(outImg);
} /* end doFinalTransform */


matrix<double> TurboRegTransform::getTransformationMatrix () {
    return getTransformationMatrix(targetPoint, sourcePoint);
}

/*********************************************************************
 Refine the landmarks.
    ********************************************************************/
void TurboRegTransform::doRegistration (
) {
    std::stack<ImageStackItem> sourceImgPyramid;
    std::stack<MaskStackItem> sourceMskPyramid;
    std::stack<ImageStackItem> targetImgPyramid;
    std::stack<MaskStackItem> targetMskPyramid;
    if (!hasSourceMask()) {
        sourceImgPyramid = sourceImg->getPyramid();
        //sourceMskPyramid = null;
        targetImgPyramid = targetImg->getPyramid();
        targetMskPyramid = targetMsk->getPyramid();
    }
    else {
        sourceImgPyramid = sourceImg->getPyramid();
        sourceMskPyramid = sourceMsk->getPyramid();
        targetImgPyramid = targetImg->getPyramid();
        targetMskPyramid = targetMsk->getPyramid();
    }
    pyramidDepth = targetImg->getPyramidDepth();
    iterationPower = (int)pow(
            (double)ITERATION_PROGRESSION, (double)pyramidDepth);

    iterationCost = 1;

    scaleBottomDownLandmarks();
    /* 1 *****************/
    while (!targetImgPyramid.empty()) {

        ImageStackItem srcImgItem = sourceImgPyramid.top();
        sourceImgPyramid.pop();

        MaskStackItem srcMskItem = sourceMskPyramid.top();
        sourceMskPyramid.pop();

        ImageStackItem tarImgItem = targetImgPyramid.top();
        targetImgPyramid.pop();

        MaskStackItem tarMskItem = targetMskPyramid.top();
        targetMskPyramid.pop();

        iterationPower /= ITERATION_PROGRESSION;

        if (transformation == BILINEAR) {
            inNx = srcImgItem.halfWidth; //((Integer)sourceImgPyramid.pop()).intValue();
            inNy = srcImgItem.halfHeight; //((Integer)sourceImgPyramid.pop()).intValue();
            inImg = srcImgItem.halfImg; //(float[])sourceImgPyramid.pop();
            if (!hasSourceMask()) {
                inMsk.clear();
            }
            else {
                inMsk = srcMskItem.halfMask; //(float[])sourceMskPyramid.pop();
            }
            outNx = tarImgItem.halfWidth; //((Integer)targetImgPyramid.pop()).intValue();
            outNy = tarImgItem.halfHeight; //((Integer)targetImgPyramid.pop()).intValue();
            outImg = tarImgItem.halfImg; //(float[])targetImgPyramid.pop();
            outMsk = tarMskItem.halfMask; //(float[])targetMskPyramid.pop();
        }
        else {
            inNx = tarImgItem.halfWidth; //((Integer)targetImgPyramid.pop()).intValue();
            inNy = tarImgItem.halfHeight; //((Integer)targetImgPyramid.pop()).intValue();
            inImg = tarImgItem.halfImg; //(float[])targetImgPyramid.pop();
            inMsk = tarMskItem.halfMask; //(float[])targetMskPyramid.pop();
            outNx = srcImgItem.halfWidth; //((Integer)sourceImgPyramid.pop()).intValue();
            outNy = srcImgItem.halfHeight; //(((Integer)sourceImgPyramid.pop()).intValue();
            outImg = srcImgItem.halfImg; //(float[])sourceImgPyramid.pop();
            xGradient = srcImgItem.xGradient; //(float[])sourceImgPyramid.pop();
            yGradient = srcImgItem.yGradient; //(float[])sourceImgPyramid.pop();
            if (!hasSourceMask()) {
                outMsk.clear();
            }
            else {
                outMsk = srcMskItem.halfMask; //(float[])sourceMskPyramid.pop();
            }
        }
        /* 2 *****************/
        twiceInNx = 2 * inNx;
        twiceInNy = 2 * inNy;
        switch (transformation) {
            case TRANSLATION: {
                targetJacobian = 1.0;
                inverseMarquardtLevenbergOptimization();
                break;
            }
            case RIGID_BODY: {
                inverseMarquardtLevenbergRigidBodyOptimization();
                break;
            }
            case SCALED_ROTATION: {
                targetJacobian = (targetPoint(0, 0) - targetPoint(1, 0))
                        * (targetPoint(0, 0) - targetPoint(1, 0))
                        + (targetPoint(0, 1) - targetPoint(1, 1))
                        * (targetPoint(0, 1) - targetPoint(1, 1));
                inverseMarquardtLevenbergOptimization();
                break;
            }
            case AFFINE: {
                targetJacobian = (targetPoint(1, 0) - targetPoint(2, 0))
                        * targetPoint(0, 1)
                        + (targetPoint(2, 0) - targetPoint(0, 0))
                        * targetPoint(1, 1)
                        + (targetPoint(0, 0) - targetPoint(1, 0))
                        * targetPoint(2, 1);
                inverseMarquardtLevenbergOptimization();
                break;
            }
            case BILINEAR: {
                marquardtLevenbergOptimization();
                break;
            }
            default:
                break;
        }
        scaleUpLandmarks();
        sourcePh->setPoints(sourcePoint);
        iterationCost *= ITERATION_PROGRESSION;
    }
    /* 3 *****************/
    iterationPower /= ITERATION_PROGRESSION;
    if (transformation == BILINEAR) {
        inNx = sourceImg->getWidth();
        inNy = sourceImg->getHeight();
        inImg = sourceImg->getCoefficient();
        if (!hasSourceMask()) {
            inMsk.clear();
        }
        else {
            inMsk = sourceMsk->getMask();
        }
        outNx = targetImg->getWidth();
        outNy = targetImg->getHeight();
        outImg = targetImg->getImage();
        outMsk = targetMsk->getMask();
    }
    else {
        inNx = targetImg->getWidth();
        inNy = targetImg->getHeight();
        inImg = targetImg->getCoefficient();
        inMsk = targetMsk->getMask();
        outNx = sourceImg->getWidth();
        outNy = sourceImg->getHeight();
        outImg = sourceImg->getImage();
        xGradient = sourceImg->getXGradient();
        yGradient = sourceImg->getYGradient();
        if (!hasSourceMask()) {
            outMsk.clear();
        }
        else {
            outMsk = sourceMsk->getMask();
        }
    }
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;

    switch (transformation) {
        case RIGID_BODY: {
            inverseMarquardtLevenbergRigidBodyOptimization();
            break;
        }
        case TRANSLATION:
        case SCALED_ROTATION:
        case AFFINE: {
            inverseMarquardtLevenbergOptimization();
            break;
        }
        case BILINEAR: {
            marquardtLevenbergOptimization();
            break;
        }
    }
    sourcePh->setPoints(sourcePoint);
    iterationPower = (int)pow(
            (double)ITERATION_PROGRESSION, (double)pyramidDepth);

} /* end doRegistration */

/*********************************************************************
 Save the current landmarks into a text file and return the path
    and name of the file. Rigid format.
    @see TurboRegDialog#loadLandmarks()
    ********************************************************************/
/*std::string TurboRegTransform::saveTransformation (
        std::string filename
) {
    inNx = sourceImg->getWidth();
    inNy = sourceImg->getHeight();
    outNx = targetImg->getWidth();
    outNy = targetImg->getHeight();
    std::string path = "";
    if (filename == null) {
        Frame f = new Frame();
        FileDialog fd = new FileDialog(f, "Save landmarks",
                FileDialog.SAVE);
        filename = "landmarks.txt";
        fd.setFile(filename);
        fd.setVisible(true);
        path = fd.getDirectory();
        filename = fd.getFile();
        if ((path == null) || (filename == null)) {
            return("");
        }
    }
    try {
        FileWriter fw = new FileWriter(path + filename);
        fw.write("Transformation\n");
        switch (transformation) {
            case TRANSLATION: {
                fw.write("TRANSLATION\n");
                break;
            }
            case RIGID_BODY: {
                fw.write("RIGID_BODY\n");
                break;
            }
            case SCALED_ROTATION: {
                fw.write("SCALED_ROTATION\n");
                break;
            }
            case AFFINE: {
                fw.write("AFFINE\n");
                break;
            }
            case BILINEAR: {
                fw.write("BILINEAR\n");
                break;
            }
        }
        fw.write("\n");
        fw.write("Source size\n");
        fw.write(inNx + "\t" + inNy + "\n");
        fw.write("\n");
        fw.write("Target size\n");
        fw.write(outNx + "\t" + outNy + "\n");
        fw.write("\n");
        fw.write("Refined source landmarks\n");
        if (transformation == RIGID_BODY) {
            for (int i = 0; (i < transformation); i++) {
                fw.write(sourcePoint(i, 0) + "\t" + sourcePoint(i, 1) + "\n");
            }
        }
        else {
            for (int i = 0; (i < (transformation / 2)); i++) {
                fw.write(sourcePoint(i, 0) + "\t" + sourcePoint(i, 1) + "\n");
            }
        }
        fw.write("\n");
        fw.write("Target landmarks\n");
        if (transformation == RIGID_BODY) {
            for (int i = 0; (i < transformation); i++) {
                fw.write(targetPoint(i, 0) + "\t" + targetPoint(i, 1) + "\n");
            }
        }
        else {
            for (int i = 0; (i < (transformation / 2)); i++) {
                fw.write(targetPoint(i, 0) + "\t" + targetPoint(i, 1) + "\n");
            }
        }
        fw.close();
    } catch (IOException e) {
        IJ.log(
                "IOException exception " + e.getMessage());
    } catch (SecurityException e) {
        IJ.log(
                "Security exception " + e.getMessage());
    }
    return(path + filename);
}*/ /* end saveTransformation */

/*....................................................................
methods
....................................................................*/
/*------------------------------------------------------------------*/
void TurboRegTransform::affineTransform (
        matrix<double> &m
) {
    double yx;
    double yy;
    double x0;
    double y0;
    int xMsk;
    int yMsk;
    int k = 0;

    yx = m(0, 0);
    yy = m(1, 0);
    for (int v = 0; (v < outNy); v++) {
        x0 = yx;
        y0 = yy;
        for (int u = 0; (u < outNx); u++) {
            x = x0;
            y = y0;
            xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= xMsk) && (xMsk < inNx) && (0 <= yMsk) && (yMsk < inNy)) {
                xMsk += yMsk * inNx;
                if (accelerated) {
                    outImg[k++] = inImg[xMsk];
                }
                else {
                    xIndexes();
                    yIndexes();
                    x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                    y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                    xWeights();
                    yWeights();
                    outImg[k++] = (float)interpolate();
                }
            }
            else {
                outImg[k++] = 0.0F;
            }
            x0 += m(0, 1);
            y0 += m(1, 1);
        }
        yx += m(0, 2);
        yy += m(1, 2);
    }
} /* affineTransform */

/*------------------------------------------------------------------*/
void TurboRegTransform::affineTransform (
        matrix<double> &m,
        std::vector<double> &outMsk
) {
    double yx;
    double yy;
    double x0;
    double y0;
    int xMsk;
    int yMsk;
    int k = 0;

    yx = m(0, 0);
    yy = m(1, 0);
    for (int v = 0; (v < outNy); v++) {
        x0 = yx;
        y0 = yy;
        for (int u = 0; (u < outNx); u++) {
            x = x0;
            y = y0;
            xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= xMsk) && (xMsk < inNx) && (0 <= yMsk) && (yMsk < inNy)) {
                xMsk += yMsk * inNx;
                if (accelerated) {
                    outImg[k] = inImg[xMsk];
                }
                else {
                    xIndexes();
                    yIndexes();
                    x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                    y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                    xWeights();
                    yWeights();
                    outImg[k] = (float)interpolate();
                }
                outMsk[k++] = inMsk[xMsk];
            }
            else {
                outImg[k] = 0.0F;
                outMsk[k++] = 0.0F;
            }
            x0 += m(0, 1);
            y0 += m(1, 1);
        }
        yx += m(0, 2);
        yy += m(1, 2);

    }
} /* affineTransform */

/*------------------------------------------------------------------*/
void TurboRegTransform::bilinearTransform (
        matrix<double> &m
) {
    double yx;
    double yy;
    double yxy;
    double yyy;
    double x0;
    double y0;
    int xMsk;
    int yMsk;
    int k = 0;

    yx = m(0, 0);
    yy = m(1, 0);
    yxy = 0.0;
    yyy = 0.0;
    for (int v = 0; (v < outNy); v++) {
        x0 = yx;
        y0 = yy;
        for (int u = 0; (u < outNx); u++) {
            x = x0;
            y = y0;
            xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= xMsk) && (xMsk < inNx) && (0 <= yMsk) && (yMsk < inNy)) {
                xMsk += yMsk * inNx;
                if (accelerated) {
                    outImg[k++] = inImg[xMsk];
                }
                else {
                    xIndexes();
                    yIndexes();
                    x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                    y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                    xWeights();
                    yWeights();
                    outImg[k++] = (float)interpolate();
                }
            }
            else {
                outImg[k++] = 0.0F;
            }
            x0 += m(0, 1) + yxy;
            y0 += m(1, 1) + yyy;
        }
        yx += m(0, 2);
        yy += m(1, 2);
        yxy += m(0, 3);
        yyy += m(1, 3);

    }

} /* bilinearTransform */

/*------------------------------------------------------------------*/
void TurboRegTransform::bilinearTransform (
        matrix<double> &m,
        std::vector<double> &outMsk
) {
    double yx;
    double yy;
    double yxy;
    double yyy;
    double x0;
    double y0;
    int xMsk;
    int yMsk;
    int k = 0;

    yx = m(0, 0);
    yy = m(1, 0);
    yxy = 0.0;
    yyy = 0.0;
    for (int v = 0; (v < outNy); v++) {
        x0 = yx;
        y0 = yy;
        for (int u = 0; (u < outNx); u++) {
            x = x0;
            y = y0;
            xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= xMsk) && (xMsk < inNx) && (0 <= yMsk) && (yMsk < inNy)) {
                xMsk += yMsk * inNx;
                if (accelerated) {
                    outImg[k] = inImg[xMsk];
                }
                else {
                    xIndexes();
                    yIndexes();
                    x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                    y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                    xWeights();
                    yWeights();
                    outImg[k] = (float)interpolate();
                }
                outMsk[k++] = inMsk[xMsk];
            }
            else {
                outImg[k] = 0.0F;
                outMsk[k++] = 0.0F;
            }
            x0 += m(0, 1) + yxy;
            y0 += m(1, 1) + yyy;
        }
        yx += m(0, 2);
        yy += m(1, 2);
        yxy += m(0, 3);
        yyy += m(1, 3);

    }

} /* bilinearTransform */

/*------------------------------------------------------------------*/
void TurboRegTransform::computeBilinearGradientConstants (
) {
    double u1 = targetPoint(0, 0);
    double u2 = targetPoint(1, 0);
    double u3 = targetPoint(2, 0);
    double u4 = targetPoint(3, 0);
    double v1 = targetPoint(0, 1);
    double v2 = targetPoint(1, 1);
    double v3 = targetPoint(2, 1);
    double v4 = targetPoint(3, 1);
    double v12 = v1 - v2;
    double v13 = v1 - v3;
    double v14 = v1 - v4;
    double v23 = v2 - v3;
    double v24 = v2 - v4;
    double v34 = v3 - v4;
    double uv12 = u1 * u2 * v12;
    double uv13 = u1 * u3 * v13;
    double uv14 = u1 * u4 * v14;
    double uv23 = u2 * u3 * v23;
    double uv24 = u2 * u4 * v24;
    double uv34 = u3 * u4 * v34;
    double det = uv12 * v34 - uv13 * v24 + uv14 * v23 + uv23 * v14
            - uv24 * v13 + uv34 * v12;
    c0 = (-uv34 * v2 + uv24 * v3 - uv23 * v4) / det;
    c0u = (u3 * v3 * v24 - u2 * v2 * v34 - u4 * v4 * v23) / det;
    c0v = (uv23 - uv24 + uv34) / det;
    c0uv = (u4 * v23 - u3 * v24 + u2 * v34) / det;
    c1 = (uv34 * v1 - uv14 * v3 + uv13 * v4) / det;
    c1u = (-u3 * v3 * v14 + u1 * v1 * v34 + u4 * v4 * v13) / det;
    c1v = (-uv13 + uv14 - uv34) / det;
    c1uv = (-u4 * v13 + u3 * v14 - u1 * v34) / det;
    c2 = (-uv24 * v1 + uv14 * v2 - uv12 * v4) / det;
    c2u = (u2 * v2 * v14 - u1 * v1 * v24 - u4 * v4 * v12) / det;
    c2v = (uv12 - uv14 + uv24) / det;
    c2uv = (u4 * v12 - u2 * v14 + u1 * v24) / det;
    c3 = (uv23 * v1 - uv13 * v2 + uv12 * v3) / det;
    c3u = (-u2 * v2 * v13 + u1 * v1 * v23 + u3 * v3 * v12) / det;
    c3v = (-uv12 + uv13 - uv23) / det;
    c3uv = (-u3 * v1 + u2 * v13 + u3 * v2 - u1 * v23) / det;
} /* end computeBilinearGradientConstants */

/*------------------------------------------------------------------*/
double TurboRegTransform::getAffineMeanSquares (
        matrix<double> &sourcePoint,
        matrix<double> &m
) {
    double u1 = sourcePoint(0, 0);
    double u2 = sourcePoint(1, 0);
    double u3 = sourcePoint(2, 0);
    double v1 = sourcePoint(0, 1);
    double v2 = sourcePoint(1, 1);
    double v3 = sourcePoint(2, 1);
    double uv32 = u3 * v2 - u2 * v3;
    double uv21 = u2 * v1 - u1 * v2;
    double uv13 = u1 * v3 - u3 * v1;
    double det = uv32 + uv21 + uv13;
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    return(meanSquares / ((double)area * std::abs(det / targetJacobian)));
} /* getAffineMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getAffineMeanSquares (
        matrix<double> &sourcePoint,
        matrix<double> &m,
        std::vector<double> &gradient
) {
    double u1 = sourcePoint(0, 0);
    double u2 = sourcePoint(1, 0);
    double u3 = sourcePoint(2, 0);
    double v1 = sourcePoint(0, 1);
    double v2 = sourcePoint(1, 1);
    double v3 = sourcePoint(2, 1);
    double uv32 = u3 * v2 - u2 * v3;
    double uv21 = u2 * v1 - u1 * v2;
    double uv13 = u1 * v3 - u3 * v1;
    double det = uv32 + uv21 + uv13;
    double u12 = (u1 - u2) /det;
    double u23 = (u2 - u3) /det;
    double u31 = (u3 - u1) /det;
    double v12 = (v1 - v2) /det;
    double v23 = (v2 - v3) /det;
    double v31 = (v3 - v1) /det;
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    double g0;
    double g1;
    double g2;
    double dx0;
    double dx1;
    double dx2;
    double dy0;
    double dy1;
    double dy2;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    uv32 /= det;
    uv21 /= det;
    uv13 /= det;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
    }
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        g0 = u23 * (double)v - v23 * (double)u + uv32;
                        g1 = u31 * (double)v - v31 * (double)u + uv13;
                        g2 = u12 * (double)v - v12 * (double)u + uv21;
                        dx0 = xGradient[k] * g0;
                        dy0 = yGradient[k] * g0;
                        dx1 = xGradient[k] * g1;
                        dy1 = yGradient[k] * g1;
                        dx2 = xGradient[k] * g2;
                        dy2 = yGradient[k] * g2;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        gradient[4] += difference * dx2;
                        gradient[5] += difference * dy2;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        g0 = u23 * (double)v - v23 * (double)u + uv32;
                        g1 = u31 * (double)v - v31 * (double)u + uv13;
                        g2 = u12 * (double)v - v12 * (double)u + uv21;
                        dx0 = xGradient[k] * g0;
                        dy0 = yGradient[k] * g0;
                        dx1 = xGradient[k] * g1;
                        dy1 = yGradient[k] * g1;
                        dx2 = xGradient[k] * g2;
                        dy2 = yGradient[k] * g2;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        gradient[4] += difference * dx2;
                        gradient[5] += difference * dy2;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    return(meanSquares / ((double)area * std::abs(det / targetJacobian)));
} /* getAffineMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getAffineMeanSquares (
        matrix<double> &sourcePoint,
        matrix<double> &m,
        matrix<double> &hessian,
        std::vector<double> &gradient
) {
    double u1 = sourcePoint(0, 0);
    double u2 = sourcePoint(1, 0);
    double u3 = sourcePoint(2, 0);
    double v1 = sourcePoint(0, 1);
    double v2 = sourcePoint(1, 1);
    double v3 = sourcePoint(2, 1);
    double uv32 = u3 * v2 - u2 * v3;
    double uv21 = u2 * v1 - u1 * v2;
    double uv13 = u1 * v3 - u3 * v1;
    double det = uv32 + uv21 + uv13;
    double u12 = (u1 - u2) /det;
    double u23 = (u2 - u3) /det;
    double u31 = (u3 - u1) /det;
    double v12 = (v1 - v2) /det;
    double v23 = (v2 - v3) /det;
    double v31 = (v3 - v1) /det;
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    double g0;
    double g1;
    double g2;
    double dx0;
    double dx1;
    double dx2;
    double dy0;
    double dy1;
    double dy2;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    uv32 /= det;
    uv21 /= det;
    uv13 /= det;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
        for (int j = 0; (j < transformation); j++) {
            hessian(i, j) = 0.0;
        }
    }
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        g0 = u23 * (double)v - v23 * (double)u + uv32;
                        g1 = u31 * (double)v - v31 * (double)u + uv13;
                        g2 = u12 * (double)v - v12 * (double)u + uv21;
                        dx0 = xGradient[k] * g0;
                        dy0 = yGradient[k] * g0;
                        dx1 = xGradient[k] * g1;
                        dy1 = yGradient[k] * g1;
                        dx2 = xGradient[k] * g2;
                        dy2 = yGradient[k] * g2;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        gradient[4] += difference * dx2;
                        gradient[5] += difference * dy2;
                        hessian(0, 0) += dx0 * dx0;
                        hessian(0, 1) += dx0 * dy0;
                        hessian(0, 2) += dx0 * dx1;
                        hessian(0, 3) += dx0 * dy1;
                        hessian(0, 4) += dx0 * dx2;
                        hessian(0, 5) += dx0 * dy2;
                        hessian(1, 1) += dy0 * dy0;
                        hessian(1, 2) += dy0 * dx1;
                        hessian(1, 3) += dy0 * dy1;
                        hessian(1, 4) += dy0 * dx2;
                        hessian(1, 5) += dy0 * dy2;
                        hessian(2, 2) += dx1 * dx1;
                        hessian(2, 3) += dx1 * dy1;
                        hessian(2, 4) += dx1 * dx2;
                        hessian(2, 5) += dx1 * dy2;
                        hessian(3, 3) += dy1 * dy1;
                        hessian(3, 4) += dy1 * dx2;
                        hessian(3, 5) += dy1 * dy2;
                        hessian(4, 4) += dx2 * dx2;
                        hessian(4, 5) += dx2 * dy2;
                        hessian(5, 5) += dy2 * dy2;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        g0 = u23 * (double)v - v23 * (double)u + uv32;
                        g1 = u31 * (double)v - v31 * (double)u + uv13;
                        g2 = u12 * (double)v - v12 * (double)u + uv21;
                        dx0 = xGradient[k] * g0;
                        dy0 = yGradient[k] * g0;
                        dx1 = xGradient[k] * g1;
                        dy1 = yGradient[k] * g1;
                        dx2 = xGradient[k] * g2;
                        dy2 = yGradient[k] * g2;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        gradient[4] += difference * dx2;
                        gradient[5] += difference * dy2;
                        hessian(0, 0) += dx0 * dx0;
                        hessian(0, 1) += dx0 * dy0;
                        hessian(0, 2) += dx0 * dx1;
                        hessian(0, 3) += dx0 * dy1;
                        hessian(0, 4) += dx0 * dx2;
                        hessian(0, 5) += dx0 * dy2;
                        hessian(1, 1) += dy0 * dy0;
                        hessian(1, 2) += dy0 * dx1;
                        hessian(1, 3) += dy0 * dy1;
                        hessian(1, 4) += dy0 * dx2;
                        hessian(1, 5) += dy0 * dy2;
                        hessian(2, 2) += dx1 * dx1;
                        hessian(2, 3) += dx1 * dy1;
                        hessian(2, 4) += dx1 * dx2;
                        hessian(2, 5) += dx1 * dy2;
                        hessian(3, 3) += dy1 * dy1;
                        hessian(3, 4) += dy1 * dx2;
                        hessian(3, 5) += dy1 * dy2;
                        hessian(4, 4) += dx2 * dx2;
                        hessian(4, 5) += dx2 * dy2;
                        hessian(5, 5) += dy2 * dy2;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian(i, j) = hessian(j, i);
        }
    }
    return(meanSquares / ((double)area * std::abs(det / targetJacobian)));
} /* getAffineMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getBilinearMeanSquares (
        matrix<double> &m
) {
    double yx;
    double yy;
    double yxy;
    double yyy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    if (inMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        yxy = 0.0;
        yyy = 0.0;
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                if (outMsk[k] != 0.0F) {
                    x = x0;
                    y = y0;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)
                            && (0 <= yMsk) && (yMsk < inNy)) {
                        xIndexes();
                        yIndexes();
                        area++;
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = interpolate() - (double)outImg[k];
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1) + yxy;
                y0 += m(1, 1) + yyy;
            }
            yx += m(0, 2);
            yy += m(1, 2);
            yxy += m(0, 3);
            yyy += m(1, 3);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        yxy = 0.0;
        yyy = 0.0;
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    xMsk += yMsk * inNx;
                    if ((outMsk[k] * inMsk[xMsk]) != 0.0F) {
                        xIndexes();
                        yIndexes();
                        area++;
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = interpolate() - (double)outImg[k];
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1) + yxy;
                y0 += m(1, 1) + yyy;
            }
            yx += m(0, 2);
            yy += m(1, 2);
            yxy += m(0, 3);
            yyy += m(1, 3);
        }
    }
    return(meanSquares / (double)area);
} /* getBilinearMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getBilinearMeanSquares (
        matrix<double> &m,
        matrix<double> &hessian,
        std::vector<double> &gradient
) {
    double yx;
    double yy;
    double yxy;
    double yyy;
    double x0;
    double y0;
    double uv;
    double xGradient;
    double yGradient;
    double difference;
    double meanSquares = 0.0;
    double g0;
    double g1;
    double g2;
    double g3;
    double dx0;
    double dx1;
    double dx2;
    double dx3;
    double dy0;
    double dy1;
    double dy2;
    double dy3;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    computeBilinearGradientConstants();
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
        for (int j = 0; (j < transformation); j++) {
            hessian(i, j) = 0.0;
        }
    }
    if (inMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        yxy = 0.0;
        yyy = 0.0;
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                if (outMsk[k] != 0.0F) {
                    x = x0;
                    y = y0;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)
                            && (0 <= yMsk) && (yMsk < inNy)) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xDxWeights();
                        yDyWeights();
                        difference = interpolate() - (double)outImg[k];
                        meanSquares += difference * difference;
                        xGradient = interpolateDx();
                        yGradient = interpolateDy();
                        uv = (double)u * (double)v;
                        g0 = c0uv * uv + c0u * (double)u + c0v * (double)v + c0;
                        g1 = c1uv * uv + c1u * (double)u + c1v * (double)v + c1;
                        g2 = c2uv * uv + c2u * (double)u + c2v * (double)v + c2;
                        g3 = c3uv * uv + c3u * (double)u + c3v * (double)v + c3;
                        dx0 = xGradient * g0;
                        dy0 = yGradient * g0;
                        dx1 = xGradient * g1;
                        dy1 = yGradient * g1;
                        dx2 = xGradient * g2;
                        dy2 = yGradient * g2;
                        dx3 = xGradient * g3;
                        dy3 = yGradient * g3;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        gradient[4] += difference * dx2;
                        gradient[5] += difference * dy2;
                        gradient[6] += difference * dx3;
                        gradient[7] += difference * dy3;
                        hessian(0, 0) += dx0 * dx0;
                        hessian(0, 1) += dx0 * dy0;
                        hessian(0, 2) += dx0 * dx1;
                        hessian(0, 3) += dx0 * dy1;
                        hessian(0, 4) += dx0 * dx2;
                        hessian(0, 5) += dx0 * dy2;
                        hessian(0, 6) += dx0 * dx3;
                        hessian(0, 7) += dx0 * dy3;
                        hessian(1, 1) += dy0 * dy0;
                        hessian(1, 2) += dy0 * dx1;
                        hessian(1, 3) += dy0 * dy1;
                        hessian(1, 4) += dy0 * dx2;
                        hessian(1, 5) += dy0 * dy2;
                        hessian(1, 6) += dy0 * dx3;
                        hessian(1, 7) += dy0 * dy3;
                        hessian(2, 2) += dx1 * dx1;
                        hessian(2, 3) += dx1 * dy1;
                        hessian(2, 4) += dx1 * dx2;
                        hessian(2, 5) += dx1 * dy2;
                        hessian(2, 6) += dx1 * dx3;
                        hessian(2, 7) += dx1 * dy3;
                        hessian(3, 3) += dy1 * dy1;
                        hessian(3, 4) += dy1 * dx2;
                        hessian(3, 5) += dy1 * dy2;
                        hessian(3, 6) += dy1 * dx3;
                        hessian(3, 7) += dy1 * dy3;
                        hessian(4, 4) += dx2 * dx2;
                        hessian(4, 5) += dx2 * dy2;
                        hessian(4, 6) += dx2 * dx3;
                        hessian(4, 7) += dx2 * dy3;
                        hessian(5, 5) += dy2 * dy2;
                        hessian(5, 6) += dy2 * dx3;
                        hessian(5, 7) += dy2 * dy3;
                        hessian(6, 6) += dx3 * dx3;
                        hessian(6, 7) += dx3 * dy3;
                        hessian(7, 7) += dy3 * dy3;
                    }
                }
                x0 += m(0, 1) + yxy;
                y0 += m(1, 1) + yyy;
            }
            yx += m(0, 2);
            yy += m(1, 2);
            yxy += m(0, 3);
            yyy += m(1, 3);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        yxy = 0.0;
        yyy = 0.0;
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    xMsk += yMsk * inNx;
                    if ((outMsk[k] * inMsk[xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xDxWeights();
                        yDyWeights();
                        difference = interpolate() - (double)outImg[k];
                        meanSquares += difference * difference;
                        xGradient = interpolateDx();
                        yGradient = interpolateDy();
                        uv = (double)u * (double)v;
                        g0 = c0uv * uv + c0u * (double)u + c0v * (double)v + c0;
                        g1 = c1uv * uv + c1u * (double)u + c1v * (double)v + c1;
                        g2 = c2uv * uv + c2u * (double)u + c2v * (double)v + c2;
                        g3 = c3uv * uv + c3u * (double)u + c3v * (double)v + c3;
                        dx0 = xGradient * g0;
                        dy0 = yGradient * g0;
                        dx1 = xGradient * g1;
                        dy1 = yGradient * g1;
                        dx2 = xGradient * g2;
                        dy2 = yGradient * g2;
                        dx3 = xGradient * g3;
                        dy3 = yGradient * g3;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        gradient[4] += difference * dx2;
                        gradient[5] += difference * dy2;
                        gradient[6] += difference * dx3;
                        gradient[7] += difference * dy3;
                        hessian(0, 0) += dx0 * dx0;
                        hessian(0, 1) += dx0 * dy0;
                        hessian(0, 2) += dx0 * dx1;
                        hessian(0, 3) += dx0 * dy1;
                        hessian(0, 4) += dx0 * dx2;
                        hessian(0, 5) += dx0 * dy2;
                        hessian(0, 6) += dx0 * dx3;
                        hessian(0, 7) += dx0 * dy3;
                        hessian(1, 1) += dy0 * dy0;
                        hessian(1, 2) += dy0 * dx1;
                        hessian(1, 3) += dy0 * dy1;
                        hessian(1, 4) += dy0 * dx2;
                        hessian(1, 5) += dy0 * dy2;
                        hessian(1, 6) += dy0 * dx3;
                        hessian(1, 7) += dy0 * dy3;
                        hessian(2, 2) += dx1 * dx1;
                        hessian(2, 3) += dx1 * dy1;
                        hessian(2, 4) += dx1 * dx2;
                        hessian(2, 5) += dx1 * dy2;
                        hessian(2, 6) += dx1 * dx3;
                        hessian(2, 7) += dx1 * dy3;
                        hessian(3, 3) += dy1 * dy1;
                        hessian(3, 4) += dy1 * dx2;
                        hessian(3, 5) += dy1 * dy2;
                        hessian(3, 6) += dy1 * dx3;
                        hessian(3, 7) += dy1 * dy3;
                        hessian(4, 4) += dx2 * dx2;
                        hessian(4, 5) += dx2 * dy2;
                        hessian(4, 6) += dx2 * dx3;
                        hessian(4, 7) += dx2 * dy3;
                        hessian(5, 5) += dy2 * dy2;
                        hessian(5, 6) += dy2 * dx3;
                        hessian(5, 7) += dy2 * dy3;
                        hessian(6, 6) += dx3 * dx3;
                        hessian(6, 7) += dx3 * dy3;
                        hessian(7, 7) += dy3 * dy3;
                    }
                }
                x0 += m(0, 1) + yxy;
                y0 += m(1, 1) + yyy;
            }
            yx += m(0, 2);
            yy += m(1, 2);
            yxy += m(0, 3);
            yyy += m(1, 3);
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian(i, j) = hessian(j, i);
        }
    }
    return(meanSquares / (double)area);
} /* getBilinearMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getRigidBodyMeanSquares (
        matrix<double> &m
) {
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    return(meanSquares / (double)area);
} /* getRigidBodyMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getRigidBodyMeanSquares (
        matrix<double> &m,
        std::vector<double> &gradient

) {
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
    }
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        gradient[0] += difference * (yGradient[k] * (double)u
                                - xGradient[k] * (double)v);
                        gradient[1] += difference * xGradient[k];
                        gradient[2] += difference * yGradient[k];
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        gradient[0] += difference * (yGradient[k] * (double)u
                                - xGradient[k] * (double)v);
                        gradient[1] += difference * xGradient[k];
                        gradient[2] += difference * yGradient[k];
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    return(meanSquares / (double)area);
} /* getRigidBodyMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getRigidBodyMeanSquares (
        matrix<double> &m,
        matrix<double> &hessian,
        std::vector<double> &gradient
) {
    double yx;
    double yy;
    double x0;
    double y0;
    double dTheta;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
        for (int j = 0; (j < transformation); j++) {
            hessian(i, j) = 0.0;
        }
    }
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        dTheta = yGradient[k] * (double)u
                                - xGradient[k] * (double)v;
                        gradient[0] += difference * dTheta;
                        gradient[1] += difference * xGradient[k];
                        gradient[2] += difference * yGradient[k];
                        hessian(0, 0) += dTheta * dTheta;
                        hessian(0, 1) += dTheta * xGradient[k];
                        hessian(0, 2) += dTheta * yGradient[k];
                        hessian(1, 1) += xGradient[k] * xGradient[k];
                        hessian(1, 2) += xGradient[k] * yGradient[k];
                        hessian(2, 2) += yGradient[k] * yGradient[k];
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        dTheta = yGradient[k] * (double)u
                                - xGradient[k] * (double)v;
                        gradient[0] += difference * dTheta;
                        gradient[1] += difference * xGradient[k];
                        gradient[2] += difference * yGradient[k];
                        hessian(0, 0) += dTheta * dTheta;
                        hessian(0, 1) += dTheta * xGradient[k];
                        hessian(0, 2) += dTheta * yGradient[k];
                        hessian(1, 1) += xGradient[k] * xGradient[k];
                        hessian(1, 2) += xGradient[k] * yGradient[k];
                        hessian(2, 2) += yGradient[k] * yGradient[k];
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian(i, j) = hessian(j, i);
        }
    }
    return(meanSquares / (double)area);
} /* getRigidBodyMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getScaledRotationMeanSquares (
        matrix<double> &sourcePoint,
        matrix<double> &m
) {
    double u1 = sourcePoint(0, 0);
    double u2 = sourcePoint(1, 0);
    double v1 = sourcePoint(0, 1);
    double v2 = sourcePoint(1, 1);
    double u12 = u1 - u2;
    double v12 = v1 - v2;
    double uv2 = u12 * u12 + v12 * v12;
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    return(meanSquares / ((double)area * uv2 / targetJacobian));
} /* getScaledRotationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getScaledRotationMeanSquares (
        matrix<double> &sourcePoint,
        matrix<double> &m,
        std::vector<double> &gradient
) {
    double u1 = sourcePoint(0, 0);
    double u2 = sourcePoint(1, 0);
    double v1 = sourcePoint(0, 1);
    double v2 = sourcePoint(1, 1);
    double u12 = u1 - u2;
    double v12 = v1 - v2;
    double uv2 = u12 * u12 + v12 * v12;
    double c = 0.5 * (u2 * v1 - u1 * v2) / uv2;
    double c1 = u12 / uv2;
    double c2 = v12 / uv2;
    double c3 = (uv2 - u12 * v12) / uv2;
    double c4 = (uv2 + u12 * v12) / uv2;
    double c5 = c + u1 * c1 + u2 * c2;
    double c6 = c * (u12 * u12 - v12 * v12) / uv2;
    double c7 = c1 * c4;
    double c8 = c1 - c2 - c1 * c2 * v12;
    double c9 = c1 + c2 - c1 * c2 * u12;
    double c0 = c2 * c3;
    double dgxx0 = c1 * u2 + c2 * v2;
    double dgyx0 = 2.0 * c;
    double dgxx1 = c5 + c6;
    double dgyy1 = c5 - c6;
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    double gxx0;
    double gxx1;
    double gxy0;
    double gxy1;
    double gyx0;
    double gyx1;
    double gyy0;
    double gyy1;
    double dx0;
    double dx1;
    double dy0;
    double dy1;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
    }
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        gxx0 = (double)u * c1 + (double)v * c2 - dgxx0;
                        gyx0 = (double)v * c1 - (double)u * c2 + dgyx0;
                        gxy0 = -gyx0;
                        gyy0 = gxx0;
                        gxx1 = (double)v * c8 - (double)u * c7 + dgxx1;
                        gyx1 = -c3 * gyx0;
                        gxy1 = c4 * gyx0;
                        gyy1 = dgyy1 - (double)u * c9 - (double)v * c0;
                        dx0 = xGradient[k] * gxx0 + yGradient[k] * gyx0;
                        dy0 = xGradient[k] * gxy0 + yGradient[k] * gyy0;
                        dx1 = xGradient[k] * gxx1 + yGradient[k] * gyx1;
                        dy1 = xGradient[k] * gxy1 + yGradient[k] * gyy1;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        gxx0 = (double)u * c1 + (double)v * c2 - dgxx0;
                        gyx0 = (double)v * c1 - (double)u * c2 + dgyx0;
                        gxy0 = -gyx0;
                        gyy0 = gxx0;
                        gxx1 = (double)v * c8 - (double)u * c7 + dgxx1;
                        gyx1 = -c3 * gyx0;
                        gxy1 = c4 * gyx0;
                        gyy1 = dgyy1 - (double)u * c9 - (double)v * c0;
                        dx0 = xGradient[k] * gxx0 + yGradient[k] * gyx0;
                        dy0 = xGradient[k] * gxy0 + yGradient[k] * gyy0;
                        dx1 = xGradient[k] * gxx1 + yGradient[k] * gyx1;
                        dy1 = xGradient[k] * gxy1 + yGradient[k] * gyy1;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    return(meanSquares / ((double)area * uv2 / targetJacobian));
} /* getScaledRotationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getScaledRotationMeanSquares (
        matrix<double> &sourcePoint,
        matrix<double> &m,
        matrix<double> &hessian,
        std::vector<double> &gradient
) {
    double u1 = sourcePoint(0, 0);
    double u2 = sourcePoint(1, 0);
    double v1 = sourcePoint(0, 1);
    double v2 = sourcePoint(1, 1);
    double u12 = u1 - u2;
    double v12 = v1 - v2;
    double uv2 = u12 * u12 + v12 * v12;
    double c = 0.5 * (u2 * v1 - u1 * v2) / uv2;
    double c1 = u12 / uv2;
    double c2 = v12 / uv2;
    double c3 = (uv2 - u12 * v12) / uv2;
    double c4 = (uv2 + u12 * v12) / uv2;
    double c5 = c + u1 * c1 + u2 * c2;
    double c6 = c * (u12 * u12 - v12 * v12) / uv2;
    double c7 = c1 * c4;
    double c8 = c1 - c2 - c1 * c2 * v12;
    double c9 = c1 + c2 - c1 * c2 * u12;
    double c0 = c2 * c3;
    double dgxx0 = c1 * u2 + c2 * v2;
    double dgyx0 = 2.0 * c;
    double dgxx1 = c5 + c6;
    double dgyy1 = c5 - c6;
    double yx;
    double yy;
    double x0;
    double y0;
    double difference;
    double meanSquares = 0.0;
    double gxx0;
    double gxx1;
    double gxy0;
    double gxy1;
    double gyx0;
    double gyx1;
    double gyy0;
    double gyy1;
    double dx0;
    double dx1;
    double dy0;
    double dy1;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
        for (int j = 0; (j < transformation); j++) {
            hessian(i, j) = 0.0;
        }
    }
    if (outMsk.empty()) {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if (inMsk[yMsk * inNx + xMsk] != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        gxx0 = (double)u * c1 + (double)v * c2 - dgxx0;
                        gyx0 = (double)v * c1 - (double)u * c2 + dgyx0;
                        gxy0 = -gyx0;
                        gyy0 = gxx0;
                        gxx1 = (double)v * c8 - (double)u * c7 + dgxx1;
                        gyx1 = -c3 * gyx0;
                        gxy1 = c4 * gyx0;
                        gyy1 = dgyy1 - (double)u * c9 - (double)v * c0;
                        dx0 = xGradient[k] * gxx0 + yGradient[k] * gyx0;
                        dy0 = xGradient[k] * gxy0 + yGradient[k] * gyy0;
                        dx1 = xGradient[k] * gxx1 + yGradient[k] * gyx1;
                        dy1 = xGradient[k] * gxy1 + yGradient[k] * gyy1;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        hessian(0, 0) += dx0 * dx0;
                        hessian(0, 1) += dx0 * dy0;
                        hessian(0, 2) += dx0 * dx1;
                        hessian(0, 3) += dx0 * dy1;
                        hessian(1, 1) += dy0 * dy0;
                        hessian(1, 2) += dy0 * dx1;
                        hessian(1, 3) += dy0 * dy1;
                        hessian(2, 2) += dx1 * dx1;
                        hessian(2, 3) += dx1 * dy1;
                        hessian(3, 3) += dy1 * dy1;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    else {
        yx = m(0, 0);
        yy = m(1, 0);
        for (int v = 0; (v < outNy); v++) {
            x0 = yx;
            y0 = yy;
            for (int u = 0; (u < outNx); u++, k++) {
                x = x0;
                y = y0;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)
                        && (0 <= yMsk) && (yMsk < inNy)) {
                    if ((outMsk[k] * inMsk[yMsk * inNx + xMsk]) != 0.0F) {
                        area++;
                        xIndexes();
                        yIndexes();
                        x -= (0.0 <= x) ? ((int)x) : ((int)x - 1);
                        y -= (0.0 <= y) ? ((int)y) : ((int)y - 1);
                        xWeights();
                        yWeights();
                        difference = (double)outImg[k] - interpolate();
                        meanSquares += difference * difference;
                        gxx0 = (double)u * c1 + (double)v * c2 - dgxx0;
                        gyx0 = (double)v * c1 - (double)u * c2 + dgyx0;
                        gxy0 = -gyx0;
                        gyy0 = gxx0;
                        gxx1 = (double)v * c8 - (double)u * c7 + dgxx1;
                        gyx1 = -c3 * gyx0;
                        gxy1 = c4 * gyx0;
                        gyy1 = dgyy1 - (double)u * c9 - (double)v * c0;
                        dx0 = xGradient[k] * gxx0 + yGradient[k] * gyx0;
                        dy0 = xGradient[k] * gxy0 + yGradient[k] * gyy0;
                        dx1 = xGradient[k] * gxx1 + yGradient[k] * gyx1;
                        dy1 = xGradient[k] * gxy1 + yGradient[k] * gyy1;
                        gradient[0] += difference * dx0;
                        gradient[1] += difference * dy0;
                        gradient[2] += difference * dx1;
                        gradient[3] += difference * dy1;
                        hessian(0, 0) += dx0 * dx0;
                        hessian(0, 1) += dx0 * dy0;
                        hessian(0, 2) += dx0 * dx1;
                        hessian(0, 3) += dx0 * dy1;
                        hessian(1, 1) += dy0 * dy0;
                        hessian(1, 2) += dy0 * dx1;
                        hessian(1, 3) += dy0 * dy1;
                        hessian(2, 2) += dx1 * dx1;
                        hessian(2, 3) += dx1 * dy1;
                        hessian(3, 3) += dy1 * dy1;
                    }
                }
                x0 += m(0, 1);
                y0 += m(1, 1);
            }
            yx += m(0, 2);
            yy += m(1, 2);
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian(i, j) = hessian(j, i);
        }
    }
    return(meanSquares / ((double)area * uv2 / targetJacobian));
} /* getScaledRotationMeanSquares */

/*------------------------------------------------------------------*/
matrix<double> TurboRegTransform::getTransformationMatrix (
        matrix<double> &fromCoord,
        matrix<double> &toCoord
) {
    matrix<double> m;
    matrix<double> a;
    std::vector<double> v;
    switch (transformation) {
        case TRANSLATION: {
            //matrix = new double[2][1];
            m.resize(2,1);
            m(0, 0) = toCoord(0, 0) - fromCoord(0, 0);
            m(1, 0) = toCoord(0, 1) - fromCoord(0, 1);
            break;
        }
        case RIGID_BODY: {
            double angle = atan2(fromCoord(2, 0) - fromCoord(1, 0),
                    fromCoord(2, 1) - fromCoord(1, 1))
                    - atan2(toCoord(2, 0) - toCoord(1, 0),
                    toCoord(2, 1) - toCoord(1, 1));
            double c = cos(angle);
            double s = sin(angle);
            //matrix = new double[2][3];
            m.resize(2,3);
            m(0, 0) = toCoord(0, 0)
                    - c * fromCoord(0, 0) + s * fromCoord(0, 1);
            m(0, 1) = c;
            m(0, 2) = -s;
            m(1, 0) = toCoord(0, 1)
                    - s * fromCoord(0, 0) - c * fromCoord(0, 1);
            m(1, 1) = s;
            m(1, 2) = c;
            break;
        }
        case SCALED_ROTATION: {
            //matrix = new double[2][3];
            m.resize(2,3);
            a.resize(3,3);
            v.resize(3);
            a(0, 0) = 1.0;
            a(0, 1) = fromCoord(0, 0);
            a(0, 2) = fromCoord(0, 1);
            a(1, 0) = 1.0;
            a(1, 1) = fromCoord(1, 0);
            a(1, 2) = fromCoord(1, 1);
            a(2, 0) = 1.0;
            a(2, 1) = fromCoord(0, 1) - fromCoord(1, 1) + fromCoord(1, 0);
            a(2, 2) = fromCoord(1, 0) + fromCoord(1, 1) - fromCoord(0, 0);
            invertGauss(a);
            v[0] = toCoord(0, 0);
            v[1] = toCoord(1, 0);
            v[2] = toCoord(0, 1) - toCoord(1, 1) + toCoord(1, 0);
            for (int i = 0; (i < 3); i++) {
                m(0, i) = 0.0;
                for (int j = 0; (j < 3); j++) {
                    m(0, i) += a(i, j) * v[j];
                }
            }
            v[0] = toCoord(0, 1);
            v[1] = toCoord(1, 1);
            v[2] = toCoord(1, 0) + toCoord(1, 1) - toCoord(0, 0);
            for (int i = 0; (i < 3); i++) {
                m(1, i) = 0.0;
                for (int j = 0; (j < 3); j++) {
                    m(1, i) += a(i, j) * v[j];
                }
            }
            break;
        }
        case AFFINE: {
            //matrix = new double[2][3];
            m.resize(2,3);
            a.resize(3,3);// = new double[3][3];
            v.resize(3);// = new double[3];
            a(0, 0) = 1.0;
            a(0, 1) = fromCoord(0, 0);
            a(0, 2) = fromCoord(0, 1);
            a(1, 0) = 1.0;
            a(1, 1) = fromCoord(1, 0);
            a(1, 2) = fromCoord(1, 1);
            a(2, 0) = 1.0;
            a(2, 1) = fromCoord(2, 0);
            a(2, 2) = fromCoord(2, 1);
            invertGauss(a);
            v[0] = toCoord(0, 0);
            v[1] = toCoord(1, 0);
            v[2] = toCoord(2, 0);
            for (int i = 0; (i < 3); i++) {
                m(0, i) = 0.0;
                for (int j = 0; (j < 3); j++) {
                    m(0, i) += a(i, j) * v[j];
                }
            }
            v[0] = toCoord(0, 1);
            v[1] = toCoord(1, 1);
            v[2] = toCoord(2, 1);
            for (int i = 0; (i < 3); i++) {
                m(1, i) = 0.0;
                for (int j = 0; (j < 3); j++) {
                    m(1, i) += a(i, j) * v[j];
                }
            }
            break;
        }
        case BILINEAR: {
            //matrix = new double[2][4];
            m.resize(2,4);
            a.resize(4,4); // = new double[4][4];
            v.resize(4); // = new double[4];
            a(0, 0) = 1.0;
            a(0, 1) = fromCoord(0, 0);
            a(0, 2) = fromCoord(0, 1);
            a(0, 3) = fromCoord(0, 0) * fromCoord(0, 1);
            a(1, 0) = 1.0;
            a(1, 1) = fromCoord(1, 0);
            a(1, 2) = fromCoord(1, 1);
            a(1, 3) = fromCoord(1, 0) * fromCoord(1, 1);
            a(2, 0) = 1.0;
            a(2, 1) = fromCoord(2, 0);
            a(2, 2) = fromCoord(2, 1);
            a(2, 3) = fromCoord(2, 0) * fromCoord(2, 1);
            a(3, 0) = 1.0;
            a(3, 1) = fromCoord(3, 0);
            a(3, 2) = fromCoord(3, 1);
            a(3, 3) = fromCoord(3, 0) * fromCoord(3, 1);
            invertGauss(a);
            v[0] = toCoord(0, 0);
            v[1] = toCoord(1, 0);
            v[2] = toCoord(2, 0);
            v[3] = toCoord(3, 0);
            for (int i = 0; (i < 4); i++) {
                m(0, i) = 0.0;
                for (int j = 0; (j < 4); j++) {
                    m(0, i) += a(i, j) * v[j];
                }
            }
            v[0] = toCoord(0, 1);
            v[1] = toCoord(1, 1);
            v[2] = toCoord(2, 1);
            v[3] = toCoord(3, 1);
            for (int i = 0; (i < 4); i++) {
                m(1, i) = 0.0;
                for (int j = 0; (j < 4); j++) {
                    m(1, i) += a(i, j) * v[j];
                }
            }
            break;
        }
    }
    return(m);
} /* end getTransformationMatrix */

/*------------------------------------------------------------------*/
double TurboRegTransform::getTranslationMeanSquares (
        matrix<double> &m
) {
    double dx = m(0, 0);
    double dy = m(1, 0);
    double dx0 = dx;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    x = dx - floor(dx);
    y = dy - floor(dy);
    xWeights();
    yWeights();
    if (outMsk.empty()) {
        for (int v = 0; (v < outNy); v++) {
            y = dy++;
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= yMsk) && (yMsk < inNy)) {
                yMsk *= inNx;
                yIndexes();
                dx = dx0;
                for (int u = 0; (u < outNx); u++, k++) {
                    x = dx++;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)) {
                        if (inMsk[yMsk + xMsk] != 0.0F) {
                            xIndexes();
                            area++;
                            difference = (double)outImg[k] - interpolate();
                            meanSquares += difference * difference;
                        }
                    }
                }
            }
            else {
                k += outNx;
            }
        }
    }
    else {
        for (int v = 0; (v < outNy); v++) {
            y = dy++;
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= yMsk) && (yMsk < inNy)) {
                yMsk *= inNx;
                yIndexes();
                dx = dx0;
                for (int u = 0; (u < outNx); u++, k++) {
                    x = dx++;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)) {
                        if ((outMsk[k] * inMsk[yMsk + xMsk]) != 0.0F) {
                            xIndexes();
                            area++;
                            difference = (double)outImg[k] - interpolate();
                            meanSquares += difference * difference;
                        }
                    }
                }
            }
            else {
                k += outNx;
            }
        }
    }
    return(meanSquares / (double)area);
} /* end getTranslationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getTranslationMeanSquares (
        matrix<double> &m,
        std::vector<double> &gradient
) {
    double dx = m(0, 0);
    double dy = m(1, 0);
    double dx0 = dx;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
    }
    x = dx - floor(dx);
    y = dy - floor(dy);
    xWeights();
    yWeights();
    if (outMsk.empty()) {
        for (int v = 0; (v < outNy); v++) {
            y = dy++;
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= yMsk) && (yMsk < inNy)) {
                yMsk *= inNx;
                yIndexes();
                dx = dx0;
                for (int u = 0; (u < outNx); u++, k++) {
                    x = dx++;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)) {
                        if (inMsk[yMsk + xMsk] != 0.0F) {
                            area++;
                            xIndexes();
                            difference = (double)outImg[k] - interpolate();
                            meanSquares += difference * difference;
                            gradient[0] += difference * xGradient[k];
                            gradient[1] += difference * yGradient[k];
                        }
                    }
                }
            }
            else {
                k += outNx;
            }
        }
    }
    else {
        for (int v = 0; (v < outNy); v++) {
            y = dy++;
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= yMsk) && (yMsk < inNy)) {
                yMsk *= inNx;
                yIndexes();
                dx = dx0;
                for (int u = 0; (u < outNx); u++, k++) {
                    x = dx++;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)) {
                        if ((outMsk[k] * inMsk[yMsk + xMsk]) != 0.0F) {
                            area++;
                            xIndexes();
                            difference = (double)outImg[k] - interpolate();
                            meanSquares += difference * difference;
                            gradient[0] += difference * xGradient[k];
                            gradient[1] += difference * yGradient[k];
                        }
                    }
                }
            }
            else {
                k += outNx;
            }
        }
    }
    return(meanSquares / (double)area);
} /* end getTranslationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getTranslationMeanSquares (
        matrix<double> &m,
        matrix<double> &hessian,
        std::vector<double> &gradient
) {
    double dx = m(0, 0);
    double dy = m(1, 0);
    double dx0 = dx;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    for (int i = 0; (i < transformation); i++) {
        gradient[i] = 0.0;
        for (int j = 0; (j < transformation); j++) {
            hessian(i, j) = 0.0;
        }
    }
    x = dx - floor(dx);
    y = dy - floor(dy);
    xWeights();
    yWeights();
    if (outMsk.empty()) {
        for (int v = 0; (v < outNy); v++) {
            y = dy++;
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= yMsk) && (yMsk < inNy)) {
                yMsk *= inNx;
                yIndexes();
                dx = dx0;
                for (int u = 0; (u < outNx); u++, k++) {
                    x = dx++;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)) {
                        if (inMsk[yMsk + xMsk] != 0.0F) {
                            area++;
                            xIndexes();
                            difference = (double)outImg[k] - interpolate();
                            meanSquares += difference * difference;
                            gradient[0] += difference * xGradient[k];
                            gradient[1] += difference * yGradient[k];
                            hessian(0, 0) += xGradient[k] * xGradient[k];
                            hessian(0, 1) += xGradient[k] * yGradient[k];
                            hessian(1, 1) += yGradient[k] * yGradient[k];
                        }
                    }
                }
            }
            else {
                k += outNx;
            }
        }
    }
    else {
        for (int v = 0; (v < outNy); v++) {
            y = dy++;
            yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
            if ((0 <= yMsk) && (yMsk < inNy)) {
                yMsk *= inNx;
                yIndexes();
                dx = dx0;
                for (int u = 0; (u < outNx); u++, k++) {
                    x = dx++;
                    xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                    if ((0 <= xMsk) && (xMsk < inNx)) {
                        if ((outMsk[k] * inMsk[yMsk + xMsk]) != 0.0F) {
                            area++;
                            xIndexes();
                            difference = (double)outImg[k] - interpolate();
                            meanSquares += difference * difference;
                            gradient[0] += difference * xGradient[k];
                            gradient[1] += difference * yGradient[k];
                            hessian(0, 0) += xGradient[k] * xGradient[k];
                            hessian(0, 1) += xGradient[k] * yGradient[k];
                            hessian(1, 1) += yGradient[k] * yGradient[k];
                        }
                    }
                }
            }
            else {
                k += outNx;
            }
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian(i, j) = hessian(j, i);
        }
    }
    return(meanSquares / (double)area);
} /* end getTranslationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::interpolate (
) {
    t = 0.0;
    for (int j = 0; (j < 4); j++) {
        s = 0.0;
        p = yIndex[j];
        for (int i = 0; (i < 4); i++) {
            s += xWeight[i] * (double)inImg[p + xIndex[i]];
        }
        t += yWeight[j] * s;
    }
    return(t);
} /* end interpolate */

/*------------------------------------------------------------------*/
double TurboRegTransform::interpolateDx (
) {
    t = 0.0;
    for (int j = 0; (j < 4); j++) {
        s = 0.0;
        p = yIndex[j];
        for (int i = 0; (i < 4); i++) {
            s += dxWeight[i] * (double)inImg[p + xIndex[i]];
        }
        t += yWeight[j] * s;
    }
    return(t);
} /* end interpolateDx */

/*------------------------------------------------------------------*/
double TurboRegTransform::interpolateDy (
) {
    t = 0.0;
    for (int j = 0; (j < 4); j++) {
        s = 0.0;
        p = yIndex[j];
        for (int i = 0; (i < 4); i++) {
            s += xWeight[i] * (double)inImg[p + xIndex[i]];
        }
        t += dyWeight[j] * s;
    }
    return(t);
} /* end interpolateDy */

/*------------------------------------------------------------------*/
void TurboRegTransform::inverseMarquardtLevenbergOptimization (
) {
    matrix<double> attempt(transformation / 2, 2); // = new double[transformation / 2][2];
    matrix<double> hessian(transformation, transformation); // = new double[transformation][transformation];
    matrix<double> pseudoHessian(transformation, transformation); // = new double[transformation][transformation];
    std::vector<double> gradient(transformation);
    matrix<double> m = getTransformationMatrix(sourcePoint, targetPoint);
    std::vector<double> update(transformation);
    double bestMeanSquares = 0.0;
    double meanSquares = 0.0;
    double lambda = FIRST_LAMBDA;
    double displacement;
    int iteration = 0;
    switch (transformation) {
        case TRANSLATION: {
            bestMeanSquares = getTranslationMeanSquares(
                    m, hessian, gradient);
            break;
        }
        case SCALED_ROTATION: {
            bestMeanSquares = getScaledRotationMeanSquares(
                    sourcePoint, m, hessian, gradient);
            break;
        }
        case AFFINE: {
            bestMeanSquares = getAffineMeanSquares(
                    sourcePoint, m, hessian, gradient);
            break;
        }
        default:
        	break;
    }
    iteration++;
    do {
        for (int k = 0; (k < transformation); k++) {
            pseudoHessian(k, k) = (1.0 + lambda) * hessian(k, k);
        }
        invertGauss(pseudoHessian);
        update = matrixMultiply(pseudoHessian, gradient);
        displacement = 0.0;
        for (int k = 0; (k < (transformation / 2)); k++) {
            attempt(k, 0) = sourcePoint(k, 0) - update[2 * k];
            attempt(k, 1) = sourcePoint(k, 1) - update[2 * k + 1];
            displacement += sqrt(update[2 * k] * update[2 * k]
                    + update[2 * k + 1] * update[2 * k + 1]);
        }
        displacement /= 0.5 * (double)transformation;
        m = getTransformationMatrix(attempt, targetPoint);
        switch (transformation) {
            case TRANSLATION: {
                if (accelerated) {
                    meanSquares = getTranslationMeanSquares(
                            m, gradient);
                }
                else {
                    meanSquares = getTranslationMeanSquares(
                            m, hessian, gradient);
                }
                break;
            }
            case SCALED_ROTATION: {
                if (accelerated) {
                    meanSquares = getScaledRotationMeanSquares(
                            attempt, m, gradient);
                }
                else {
                    meanSquares = getScaledRotationMeanSquares(
                            attempt, m, hessian, gradient);
                }
                break;
            }
            case AFFINE: {
                if (accelerated) {
                    meanSquares = getAffineMeanSquares(
                            attempt, m, gradient);
                }
                else {
                    meanSquares = getAffineMeanSquares(
                            attempt, m, hessian, gradient);
                }
                break;
            }
            default:
                break;
        }
        iteration++;
        if (meanSquares < bestMeanSquares) {
            bestMeanSquares = meanSquares;
            for (int k = 0; (k < (transformation / 2)); k++) {
                sourcePoint(k, 0) = attempt(k, 0);
                sourcePoint(k, 1) = attempt(k, 1);
            }
            lambda /= LAMBDA_MAGSTEP;
        }
        else {
            lambda *= LAMBDA_MAGSTEP;
        }

    } while ((iteration < (maxIterations * iterationPower - 1))
            && (pixelPrecision <= displacement));
    invertGauss(hessian);
    update = matrixMultiply(hessian, gradient);
    for (int k = 0; (k < (transformation / 2)); k++) {
        attempt(k, 0) = sourcePoint(k, 0) - update[2 * k];
        attempt(k, 1) = sourcePoint(k, 1) - update[2 * k + 1];
    }
    m = getTransformationMatrix(attempt, targetPoint);
    switch (transformation) {
        case TRANSLATION: {
            meanSquares = getTranslationMeanSquares(m);
            break;
        }
        case SCALED_ROTATION: {
            meanSquares = getScaledRotationMeanSquares(attempt, m);
            break;
        }
        case AFFINE: {
            meanSquares = getAffineMeanSquares(attempt, m);
            break;
        }
        default:
            break;
    }
    iteration++;
    if (meanSquares < bestMeanSquares) {
        for (int k = 0; (k < (transformation / 2)); k++) {
            sourcePoint(k, 0) = attempt(k, 0);
            sourcePoint(k, 1) = attempt(k, 1);
        }
    }

} /* end inverseMarquardtLevenbergOptimization */

/*------------------------------------------------------------------*/
void TurboRegTransform::inverseMarquardtLevenbergRigidBodyOptimization () {
    matrix<double> attempt(2, 3); //  = new double[2][3];
    matrix<double> hessian(transformation, transformation); //  = new double[transformation][transformation];
    matrix<double> pseudoHessian(transformation, transformation); //  = new double[transformation][transformation];
    std::vector<double> gradient(transformation);
    matrix<double> m = getTransformationMatrix(targetPoint, sourcePoint);
    std::vector<double> update(transformation);
    double bestMeanSquares = 0.0;
    double meanSquares = 0.0;
    double lambda = FIRST_LAMBDA;
    double angle;
    double c;
    double s;
    double displacement;
    int iteration = 0;
    for (int k = 0; (k < transformation); k++) {
        sourcePoint(k, 0) = m(0, 0) + targetPoint(k, 0) * m(0, 1)
                + targetPoint(k, 1) * m(0, 2);
        sourcePoint(k, 1) = m(1, 0) + targetPoint(k, 0) * m(1, 1)
                + targetPoint(k, 1) * m(1, 2);
    }
    m = getTransformationMatrix(sourcePoint, targetPoint);
    bestMeanSquares = getRigidBodyMeanSquares(m, hessian, gradient);
    iteration++;
    do {
        for (int k = 0; (k < transformation); k++) {
            pseudoHessian(k, k) = (1.0 + lambda) * hessian(k, k);
        }
        invertGauss(pseudoHessian);
        update = matrixMultiply(pseudoHessian, gradient);
        angle = atan2(m(0, 2), m(0, 1)) - update[0];
        attempt(0, 1) = cos(angle);
        attempt(0, 2) = sin(angle);
        attempt(1, 1) = -attempt(0, 2);
        attempt(1, 2) = attempt(0, 1);
        c = cos(update[0]);
        s = sin(update[0]);
        attempt(0, 0) = (m(0, 0) + update[1]) * c
                - (m(1, 0) + update[2]) * s;
        attempt(1, 0) = (m(0, 0) + update[1]) * s
                + (m(1, 0) + update[2]) * c;
        displacement = sqrt(update[1] * update[1] + update[2] * update[2])
                + 0.25 * sqrt((double)(inNx * inNx) + (double)(inNy * inNy))
                * std::abs(update[0]);
        if (accelerated) {
            meanSquares = getRigidBodyMeanSquares(attempt, gradient);
        }
        else {
            meanSquares = getRigidBodyMeanSquares(attempt, hessian, gradient);
        }
        iteration++;
        if (meanSquares < bestMeanSquares) {
            bestMeanSquares = meanSquares;
            for (int i = 0; (i < 2); i++) {
                for (int j = 0; (j < 3); j++) {
                    m(i, j) = attempt(i, j);
                }
            }
            lambda /= LAMBDA_MAGSTEP;
        }
        else {
            lambda *= LAMBDA_MAGSTEP;
        }

    } while ((iteration < (maxIterations * iterationPower - 1))
            && (pixelPrecision <= displacement));
    invertGauss(hessian);
    update = matrixMultiply(hessian, gradient);
    angle = atan2(m(0, 2), m(0, 1)) - update[0];
    attempt(0, 1) = cos(angle);
    attempt(0, 2) = sin(angle);
    attempt(1, 1) = -attempt(0, 2);
    attempt(1, 2) = attempt(0, 1);
    c = cos(update[0]);
    s = sin(update[0]);
    attempt(0, 0) = (m(0, 0) + update[1]) * c
            - (m(1, 0) + update[2]) * s;
    attempt(1, 0) = (m(0, 0) + update[1]) * s
            + (m(1, 0) + update[2]) * c;
    meanSquares = getRigidBodyMeanSquares(attempt);
    iteration++;
    if (meanSquares < bestMeanSquares) {
        for (int i = 0; (i < 2); i++) {
            for (int j = 0; (j < 3); j++) {
                m(i, j) = attempt(i, j);
            }
        }
    }
    for (int k = 0; (k < transformation); k++) {
        sourcePoint(k, 0) = (targetPoint(k, 0) - m(0, 0)) * m(0, 1)
                + (targetPoint(k, 1) - m(1, 0)) * m(1, 1);
        sourcePoint(k, 1) = (targetPoint(k, 0) - m(0, 0)) * m(0, 2)
                + (targetPoint(k, 1) - m(1, 0)) * m(1, 2);
    }

} /* end inverseMarquardtLevenbergRigidBodyOptimization */

/*------------------------------------------------------------------*/
void TurboRegTransform::invertGauss (
        matrix<double> &m
) {
    int n = m.nrows();
    matrix<double> inverse(n, n);
    for (int i = 0; (i < n); i++) {
        double max = m(i, 0);
        double absMax = std::abs(max);
        for (int j = 0; (j < n); j++) {
            inverse(i, j) = 0.0;
            if (absMax < std::abs(m(i, j))) {
                max = m(i, j);
                absMax = std::abs(max);
            }
        }
        inverse(i, i) = 1.0 / max;
        for (int j = 0; (j < n); j++) {
            m(i, j) /= max;
        }
    }
    for (int j = 0; (j < n); j++) {
        double max = m(j, j);
        double absMax = std::abs(max);
        int k = j;
        for (int i = j + 1; (i < n); i++) {
            if (absMax < std::abs(m(i, j))) {
                max = m(i, j);
                absMax = std::abs(max);
                k = i;
            }
        }
        if (k != j) {
            //std::vector<double> partialLine(n - j);
            //std::vector<double> fullLine(n);

            //System.arraycopy(m[j], j, partialLine, 0, n - j);
        	// copy m[k][j->n] to m[j][j->n]
            std::vector<double> partialLine(m.getPtr(j, j), m.getPtr(j, n));

            //System.arraycopy(m[k], j, m[j], j, n - j);
            // copy m[k][j->n] to m[j][j->n]
            std::copy(m.getPtr(k,j), m.getPtr(k,n), m.getPtr(j,j));

            //System.arraycopy(partialLine, 0, m[k], j, n - j);
            // copy partialLine[0->n] to m[k][j->n]
#ifdef _MSC_VER
			std::copy(partialLine.begin(), partialLine.end(), m.getPtr(k, j));
#else
			std::copy(std::begin(partialLine), std::end(partialLine), m.getPtr(k, j));
#endif


            //System.arraycopy(inverse[j], 0, fullLine, 0, n);
            // copy inverse[j][0->n] to fullLine
            std::vector<double> fullLine(inverse.getRowPtr(j), inverse.getRowPtr(j+1));

            //System.arraycopy(inverse[k], 0, inverse[j], 0, n);
            // copy inverse[k][0->n] to inverse[j][0->n]
            std::copy(inverse.getRowPtr(k), inverse.getRowPtr(k+1), inverse.getRowPtr(j));

            //System.arraycopy(fullLine, 0, inverse[k], 0, n);
            // copy fullLine to inverse[k][0->n]
#ifdef _MSC_VER
			std::copy(fullLine.begin(), fullLine.end(), inverse.getRowPtr(k));
#else
			std::copy(std::begin(fullLine), std::end(fullLine), inverse.getRowPtr(k));
#endif


        }
        for (k = 0; (k <= j); k++) {
            inverse(j, k) /= max;
        }
        for (k = j + 1; (k < n); k++) {
            m(j, k) /= max;
            inverse(j, k) /= max;
        }
        for (int i = j + 1; (i < n); i++) {
            for (k = 0; (k <= j); k++) {
                inverse(i, k) -= m(i, j) * inverse(j, k);
            }
            for (k = j + 1; (k < n); k++) {
                m(i, k) -= m(i, j) * m(j, k);
                inverse(i, k) -= m(i, j) * inverse(j, k);
            }
        }
    }
    for (int j = n - 1; (1 <= j); j--) {
        for (int i = j - 1; (0 <= i); i--) {
            for (int k = 0; (k <= j); k++) {
                inverse(i, k) -= m(i, j) * inverse(j, k);
            }
            for (int k = j + 1; (k < n); k++) {
                m(i, k) -= m(i, j) * m(j, k);
                inverse(i, k) -= m(i, j) * inverse(j, k);
            }
        }
    }
    for (int i = 0; (i < n); i++) {
        //System.arraycopy(inverse[i], 0, m[i], 0, n);
        // copy inverse[i][0->n] to m[i][0->n]
        std::copy(inverse.getRowPtr(i), inverse.getRowPtr(i+1), m.getRowPtr(i));
    }
} /* end invertGauss */

/*------------------------------------------------------------------*/
void TurboRegTransform::marquardtLevenbergOptimization (
) {
    matrix<double> attempt(transformation / 2, 2); // = new double[transformation / 2][2];
    matrix<double> hessian(transformation, transformation); // = new double[transformation][transformation];
    matrix<double> pseudoHessian(transformation, transformation); // = new double[transformation][transformation];
    std::vector<double> gradient(transformation);
    matrix<double> m = getTransformationMatrix(targetPoint, sourcePoint);
    std::vector<double> update(transformation);
    double bestMeanSquares = 0.0;
    double meanSquares = 0.0;
    double lambda = FIRST_LAMBDA;
    double displacement;
    int iteration = 0;
    bestMeanSquares = getBilinearMeanSquares(m, hessian, gradient);
    iteration++;
    do {
        for (int k = 0; (k < transformation); k++) {
            pseudoHessian(k, k) = (1.0 + lambda) * hessian(k, k);
        }
        invertGauss(pseudoHessian);
        update = matrixMultiply(pseudoHessian, gradient);
        displacement = 0.0;
        for (int k = 0; (k < (transformation / 2)); k++) {
            attempt(k, 0) = sourcePoint(k, 0) - update[2 * k];
            attempt(k, 1) = sourcePoint(k, 1) - update[2 * k + 1];
            displacement += sqrt(update[2 * k] * update[2 * k]
                    + update[2 * k + 1] * update[2 * k + 1]);
        }
        displacement /= 0.5 * (double)transformation;
        m = getTransformationMatrix(targetPoint, attempt);
        meanSquares = getBilinearMeanSquares(m, hessian, gradient);
        iteration++;
        if (meanSquares < bestMeanSquares) {
            bestMeanSquares = meanSquares;
            for (int k = 0; (k < (transformation / 2)); k++) {
                sourcePoint(k, 0) = attempt(k, 0);
                sourcePoint(k, 1) = attempt(k, 1);
            }
            lambda /= LAMBDA_MAGSTEP;
        }
        else {
            lambda *= LAMBDA_MAGSTEP;
        }

    } while ((iteration < (maxIterations * iterationPower - 1))
            && (pixelPrecision <= displacement));
    invertGauss(hessian);
    update = matrixMultiply(hessian, gradient);
    for (int k = 0; (k < (transformation / 2)); k++) {
        attempt(k, 0) = sourcePoint(k, 0) - update[2 * k];
        attempt(k, 1) = sourcePoint(k, 1) - update[2 * k + 1];
    }
    m = getTransformationMatrix(targetPoint, attempt);
    meanSquares = getBilinearMeanSquares(m);
    iteration++;
    if (meanSquares < bestMeanSquares) {
        for (int k = 0; (k < (transformation / 2)); k++) {
            sourcePoint(k, 0) = attempt(k, 0);
            sourcePoint(k, 1) = attempt(k, 1);
        }
    }

} /* end marquardtLevenbergOptimization */

/*------------------------------------------------------------------*/
std::vector<double> TurboRegTransform::matrixMultiply (
        matrix<double> &m,
        std::vector<double> &vector
) {
    std::vector<double> result(m.nrows());
    for (unsigned int i = 0; (i < m.nrows()); i++) {
        result[i] = 0.0;
        for (unsigned int j = 0; (j < vector.size()); j++) {
            result[i] += m(i, j) * vector[j];
        }
    }
    return(result);
} /* end matrixMultiply */

/*------------------------------------------------------------------*/
void TurboRegTransform::scaleBottomDownLandmarks (
) {
    for (int depth = 1; (depth < pyramidDepth); depth++) {
        if (transformation == RIGID_BODY) {
            for (int n = 0; (n < transformation); n++) {
                sourcePoint(n, 0) *= 0.5;
                sourcePoint(n, 1) *= 0.5;
                targetPoint(n, 0) *= 0.5;
                targetPoint(n, 1) *= 0.5;
            }
        }
        else {
            for (int n = 0; (n < (transformation / 2)); n++) {
                sourcePoint(n, 0) *= 0.5;
                sourcePoint(n, 1) *= 0.5;
                targetPoint(n, 0) *= 0.5;
                targetPoint(n, 1) *= 0.5;
            }
        }
    }
} /* end scaleBottomDownLandmarks */

/*------------------------------------------------------------------*/
void TurboRegTransform::scaleUpLandmarks (
) {
    if (transformation == RIGID_BODY) {
        for (int n = 0; (n < transformation); n++) {
            sourcePoint(n, 0) *= 2.0;
            sourcePoint(n, 1) *= 2.0;
            targetPoint(n, 0) *= 2.0;
            targetPoint(n, 1) *= 2.0;
        }
    }
    else {
        for (int n = 0; (n < (transformation / 2)); n++) {
            sourcePoint(n, 0) *= 2.0;
            sourcePoint(n, 1) *= 2.0;
            targetPoint(n, 0) *= 2.0;
            targetPoint(n, 1) *= 2.0;
        }
    }
} /* end scaleUpLandmarks */

/*------------------------------------------------------------------*/
void TurboRegTransform::translationTransform (
        matrix<double> &m
) {
    double dx = m(0, 0);
    double dy = m(1, 0);
    double dx0 = dx;
    int xMsk;
    int yMsk;
    x = dx - floor(dx);
    y = dy - floor(dy);
    if (!accelerated) {
        xWeights();
        yWeights();
    }
    int k = 0;

    for (int v = 0; (v < outNy); v++) {
        y = dy++;
        yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
        if ((0 <= yMsk) && (yMsk < inNy)) {
            yMsk *= inNx;
            if (!accelerated) {
                yIndexes();
            }
            dx = dx0;
            for (int u = 0; (u < outNx); u++) {
                x = dx++;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)) {
                    xMsk += yMsk;
                    if (accelerated) {
                        outImg[k++] = inImg[xMsk];
                    }
                    else {
                        xIndexes();
                        outImg[k++] = (float)interpolate();
                    }
                }
                else {
                    outImg[k++] = 0.0F;
                }
            }
        }
        else {
            for (int u = 0; (u < outNx); u++) {
                outImg[k++] = 0.0F;
            }
        }

    }

} /* translationTransform */

/*------------------------------------------------------------------*/
void TurboRegTransform::translationTransform (
        matrix<double> &m,
        std::vector<double> &outMsk
) {
    double dx = m(0, 0);
    double dy = m(1, 0);
    double dx0 = dx;
    int xMsk;
    int yMsk;
    x = dx - floor(dx);
    y = dy - floor(dy);
    if (!accelerated) {
        xWeights();
        yWeights();
    }
    int k = 0;

    for (int v = 0; (v < outNy); v++) {
        y = dy++;
        yMsk = (0.0 <= y) ? ((int)(y + 0.5)) : ((int)(y - 0.5));
        if ((0 <= yMsk) && (yMsk < inNy)) {
            yMsk *= inNx;
            if (!accelerated) {
                yIndexes();
            }
            dx = dx0;
            for (int u = 0; (u < outNx); u++, k++) {
                x = dx++;
                xMsk = (0.0 <= x) ? ((int)(x + 0.5)) : ((int)(x - 0.5));
                if ((0 <= xMsk) && (xMsk < inNx)) {
                    xMsk += yMsk;
                    if (accelerated) {
                        outImg[k] = inImg[xMsk];
                    }
                    else {
                        xIndexes();
                        outImg[k] = (float)interpolate();
                    }
                    outMsk[k] = inMsk[xMsk];
                }
                else {
                    outImg[k] = 0.0F;
                    outMsk[k] = 0.0F;
                }
            }
        }
        else {
            for (int u = 0; (u < outNx); u++, k++) {
                outImg[k] = 0.0F;
                outMsk[k] = 0.0F;
            }
        }

    }

} /* translationTransform */

/*------------------------------------------------------------------*/
void TurboRegTransform::xDxWeights (
) {
    s = 1.0 - x;
    dxWeight[0] = 0.5 * x * x;
    xWeight[0] = x * dxWeight[0] / 3.0;
    dxWeight[3] = -0.5 * s * s;
    xWeight[3] = s * dxWeight[3] / -3.0;
    dxWeight[1] = 1.0 - 2.0 * dxWeight[0] + dxWeight[3];
    xWeight[1] = 2.0 / 3.0 + (1.0 + x) * dxWeight[3];
    dxWeight[2] = 1.5 * x * (x - 4.0/ 3.0);
    xWeight[2] = 2.0 / 3.0 - (2.0 - x) * dxWeight[0];
} /* xDxWeights */

/*------------------------------------------------------------------*/
void TurboRegTransform::xIndexes (
) {
    p = (0.0 <= x) ? ((int)x + 2) : ((int)x + 1);
    for (int k = 0; (k < 4); p--, k++) {
        q = (p < 0) ? (-1 - p) : (p);
        if (twiceInNx <= q) {
            q -= twiceInNx * (q / twiceInNx);
        }
        xIndex[k] = (inNx <= q) ? (twiceInNx - 1 - q) : (q);
    }
} /* xIndexes */

/*------------------------------------------------------------------*/
void TurboRegTransform::xWeights (
) {
    s = 1.0 - x;
    xWeight[3] = s * s * s / 6.0;
    s = x * x;
    xWeight[2] = 2.0 / 3.0 - 0.5 * s * (2.0 - x);
    xWeight[0] = s * x / 6.0;
    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
} /* xWeights */

/*------------------------------------------------------------------*/
void TurboRegTransform::yDyWeights (
) {
    t = 1.0 - y;
    dyWeight[0] = 0.5 * y * y;
    yWeight[0] = y * dyWeight[0] / 3.0;
    dyWeight[3] = -0.5 * t * t;
    yWeight[3] = t * dyWeight[3] / -3.0;
    dyWeight[1] = 1.0 - 2.0 * dyWeight[0] + dyWeight[3];
    yWeight[1] = 2.0 / 3.0 + (1.0 + y) * dyWeight[3];
    dyWeight[2] = 1.5 * y * (y - 4.0/ 3.0);
    yWeight[2] = 2.0 / 3.0 - (2.0 - y) * dyWeight[0];
} /* yDyWeights */

/*------------------------------------------------------------------*/
void TurboRegTransform::yIndexes (
) {
    p = (0.0 <= y) ? ((int)y + 2) : ((int)y + 1);
    for (int k = 0; (k < 4); p--, k++) {
        q = (p < 0) ? (-1 - p) : (p);
        if (twiceInNy <= q) {
            q -= twiceInNy * (q / twiceInNy);
        }
        yIndex[k] = (inNy <= q) ? ((twiceInNy - 1 - q) * inNx) : (q * inNx);
    }
} /* yIndexes */

/*------------------------------------------------------------------*/
void TurboRegTransform::yWeights (
) {
    t = 1.0 - y;
    yWeight[3] = t * t * t / 6.0;
    t = y * y;
    yWeight[2] = 2.0 / 3.0 - 0.5 * t * (2.0 - y);
    yWeight[0] = t * y / 6.0;
    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
} /* yWeights */

void TurboRegTransform::printMatrix(matrix<double> &m) {
    for (unsigned int i = 0; i < m.nrows(); i++) {
        for (unsigned int j = 0; j < m.ncols(); j++) {
            printf("%.2f\t", m(i,j));
        }
        printf("\n");
    }
}

void TurboRegTransform::printPoints() {
    printf("target:\n");
    printMatrix(targetPoint);
    printf("source:\n");
    printMatrix(sourcePoint);
}
