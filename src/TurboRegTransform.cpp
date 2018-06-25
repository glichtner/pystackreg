#include "TurboRegTransform.h"


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
        TurboRegImage sourceImg,
        TurboRegMask sourceMsk,
        TurboRegPointHandler sourcePh,
        TurboRegImage targetImg,
        TurboRegMask targetMsk,
        TurboRegPointHandler targetPh,
        int transformation,
        bool accelerated,
        bool interactive
) {
    this->sourceImg = sourceImg;
    this->sourceMsk = sourceMsk;
    this->sourcePh = sourcePh;
    this->targetImg = targetImg;
    this->targetMsk = targetMsk;
    this->transformation = transformation;
    this->accelerated = accelerated;
    this->interactive = interactive;
    sourcePoint = sourcePh.getPoints();
    targetPoint = targetPh.getPoints();

    if (accelerated) {
        pixelPrecision = PIXEL_LOW_PRECISION;
        maxIterations = FEW_ITERATIONS;
    }
    else {
        pixelPrecision = PIXEL_HIGH_PRECISION;
        maxIterations = MANY_ITERATIONS;
    }
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
TurboRegTransform::appendTransformation (
        String pathAndFilename
) {
    outNx = targetImg.getWidth();
    outNy = targetImg.getHeight();
    inNx = sourceImg.getWidth();
    inNy = sourceImg.getHeight();
    if (pathAndFilename == null) {
        return;
    }
   /* FileWriter fw = new FileWriter(pathAndFilename, true);
    fw.write("\n");
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            fw.write("TRANSLATION\n");
            break;
        }
        case TurboRegDialog.RIGID_BODY: {
            fw.write("RIGID_BODY\n");
            break;
        }
        case TurboRegDialog.SCALED_ROTATION: {
            fw.write("SCALED_ROTATION\n");
            break;
        }
        case TurboRegDialog.AFFINE: {
            fw.write("AFFINE\n");
            break;
        }
        case TurboRegDialog.BILINEAR: {
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
    if (transformation == TurboRegDialog.RIGID_BODY) {
        for (int i = 0; (i < transformation); i++) {
            fw.write(sourcePoint[i][0] + "\t" + sourcePoint[i][1] + "\n");
        }
    }
    else {
        for (int i = 0; (i < (transformation / 2)); i++) {
            fw.write(sourcePoint[i][0] + "\t" + sourcePoint[i][1] + "\n");
        }
    }
    fw.write("\n");
    fw.write("Target landmarks\n");
    if (transformation == TurboRegDialog.RIGID_BODY) {
        for (int i = 0; (i < transformation); i++) {
            fw.write(targetPoint[i][0] + "\t" + targetPoint[i][1] + "\n");
        }
    }
    else {
        for (int i = 0; (i < (transformation / 2)); i++) {
            fw.write(targetPoint[i][0] + "\t" + targetPoint[i][1] + "\n");
        }
    }
    fw.close();*/

} /* end appendTransformation */

/*********************************************************************
 Compute the image.
    ********************************************************************/
TurboRegTransform::doBatchFinalTransform (
        float[] pixels
) {
    if (accelerated) {
        inImg = sourceImg.getImage();
    }
    else {
        inImg = sourceImg.getCoefficient();
    }
    inNx = sourceImg.getWidth();
    inNy = sourceImg.getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;
    outImg = pixels;
    outNx = targetImg.getWidth();
    outNy = targetImg.getHeight();
    double[][] matrix = getTransformationMatrix(targetPoint, sourcePoint);
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            translationTransform(matrix);
            break;
        }
        case TurboRegDialog.RIGID_BODY:
        case TurboRegDialog.SCALED_ROTATION:
        case TurboRegDialog.AFFINE: {
            affineTransform(matrix);
            break;
        }
        case TurboRegDialog.BILINEAR: {
            bilinearTransform(matrix);
            break;
        }
    }
} /* end doBatchFinalTransform */

/*********************************************************************
 Compute the image.
    ********************************************************************/
ImagePlus TurboRegTransform::doFinalTransform (
        int width,
        int height
) {
    if (accelerated) {
        inImg = sourceImg.getImage();
    }
    else {
        inImg = sourceImg.getCoefficient();
    }
    inMsk = sourceMsk.getMask();
    inNx = sourceImg.getWidth();
    inNy = sourceImg.getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;
    ImageStack is = new ImageStack(width, height);
    FloatProcessor dataFp = new FloatProcessor(width, height);
    is.addSlice("Data", dataFp);
    FloatProcessor maskFp = new FloatProcessor(width, height);
    is.addSlice("Mask", maskFp);
    ImagePlus imp = new ImagePlus("Output", is);
    imp.setSlice(1);
    outImg = (float[])dataFp.getPixels();
    imp.setSlice(2);
    float[] outMsk = (float[])maskFp.getPixels();
    outNx = imp.getWidth();
    outNy = imp.getHeight();
    double[][] matrix = getTransformationMatrix(targetPoint, sourcePoint);
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            translationTransform(matrix, outMsk);
            break;
        }
        case TurboRegDialog.RIGID_BODY:
        case TurboRegDialog.SCALED_ROTATION:
        case TurboRegDialog.AFFINE: {
            affineTransform(matrix, outMsk);
            break;
        }
        case TurboRegDialog.BILINEAR: {
            bilinearTransform(matrix, outMsk);
            break;
        }
    }
    imp.setSlice(1);
    imp.getProcessor().resetMinAndMax();
    if (interactive) {
        imp.show();
        imp.updateAndDraw();
    }
    return(imp);
} /* end doFinalTransform */

/*********************************************************************
 Compute the image.
    ********************************************************************/
public float[] doFinalTransform (
        TurboRegImage sourceImg,
        TurboRegPointHandler sourcePh,
        TurboRegImage targetImg,
        TurboRegPointHandler targetPh,
        int transformation,
        bool accelerated
) {
    this->sourceImg = sourceImg;
    this->targetImg = targetImg;
    this->sourcePh = sourcePh;
    this->transformation = transformation;
    this->accelerated = accelerated;
    sourcePoint = sourcePh.getPoints();
    targetPoint = targetPh.getPoints();
    if (accelerated) {
        inImg = sourceImg.getImage();
    }
    else {
        inImg = sourceImg.getCoefficient();
    }
    inNx = sourceImg.getWidth();
    inNy = sourceImg.getHeight();
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;
    outNx = targetImg.getWidth();
    outNy = targetImg.getHeight();
    outImg = new float[outNx * outNy];
    double[][] matrix = getTransformationMatrix(targetPoint, sourcePoint);
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            translationTransform(matrix);
            break;
        }
        case TurboRegDialog.RIGID_BODY:
        case TurboRegDialog.SCALED_ROTATION:
        case TurboRegDialog.AFFINE: {
            affineTransform(matrix);
            break;
        }
        case TurboRegDialog.BILINEAR: {
            bilinearTransform(matrix);
            break;
        }
    }
    return(outImg);
} /* end doFinalTransform */

/*********************************************************************
 Refine the landmarks.
    ********************************************************************/
TurboRegTransform::doRegistration (
) {
    Stack<?> sourceImgPyramid;
    Stack<?> sourceMskPyramid;
    Stack<?> targetImgPyramid;
    Stack<?> targetMskPyramid;
    if (sourceMsk == null) {
        sourceImgPyramid = sourceImg.getPyramid();
        sourceMskPyramid = null;
        targetImgPyramid = (Stack<?>)targetImg.getPyramid().clone();
        targetMskPyramid = (Stack<?>)targetMsk.getPyramid().clone();
    }
    else {
        sourceImgPyramid = sourceImg.getPyramid();
        sourceMskPyramid = sourceMsk.getPyramid();
        targetImgPyramid = targetImg.getPyramid();
        targetMskPyramid = targetMsk.getPyramid();
    }
    pyramidDepth = targetImg.getPyramidDepth();
    iterationPower = (int)Math.pow(
            (double)ITERATION_PROGRESSION, (double)pyramidDepth);
    TurboRegProgressBar.addWorkload(
            pyramidDepth * maxIterations * iterationPower
                    / ITERATION_PROGRESSION
                    - (iterationPower - 1) / (ITERATION_PROGRESSION - 1));
    iterationCost = 1;
    scaleBottomDownLandmarks();
    while (!targetImgPyramid.isEmpty()) {
        iterationPower /= ITERATION_PROGRESSION;
        if (transformation == TurboRegDialog.BILINEAR) {
            inNx = ((Integer)sourceImgPyramid.pop()).intValue();
            inNy = ((Integer)sourceImgPyramid.pop()).intValue();
            inImg = (float[])sourceImgPyramid.pop();
            if (sourceMskPyramid == null) {
                inMsk = null;
            }
            else {
                inMsk = (float[])sourceMskPyramid.pop();
            }
            outNx = ((Integer)targetImgPyramid.pop()).intValue();
            outNy = ((Integer)targetImgPyramid.pop()).intValue();
            outImg = (float[])targetImgPyramid.pop();
            outMsk = (float[])targetMskPyramid.pop();
        }
        else {
            inNx = ((Integer)targetImgPyramid.pop()).intValue();
            inNy = ((Integer)targetImgPyramid.pop()).intValue();
            inImg = (float[])targetImgPyramid.pop();
            inMsk = (float[])targetMskPyramid.pop();
            outNx = ((Integer)sourceImgPyramid.pop()).intValue();
            outNy = ((Integer)sourceImgPyramid.pop()).intValue();
            outImg = (float[])sourceImgPyramid.pop();
            xGradient = (float[])sourceImgPyramid.pop();
            yGradient = (float[])sourceImgPyramid.pop();
            if (sourceMskPyramid == null) {
                outMsk = null;
            }
            else {
                outMsk = (float[])sourceMskPyramid.pop();
            }
        }
        twiceInNx = 2 * inNx;
        twiceInNy = 2 * inNy;
        switch (transformation) {
            case TurboRegDialog.TRANSLATION: {
                targetJacobian = 1.0;
                inverseMarquardtLevenbergOptimization(
                        iterationPower * maxIterations - 1);
                break;
            }
            case TurboRegDialog.RIGID_BODY: {
                inverseMarquardtLevenbergRigidBodyOptimization(
                        iterationPower * maxIterations - 1);
                break;
            }
            case TurboRegDialog.SCALED_ROTATION: {
                targetJacobian = (targetPoint[0][0] - targetPoint[1][0])
                        * (targetPoint[0][0] - targetPoint[1][0])
                        + (targetPoint[0][1] - targetPoint[1][1])
                        * (targetPoint[0][1] - targetPoint[1][1]);
                inverseMarquardtLevenbergOptimization(
                        iterationPower * maxIterations - 1);
                break;
            }
            case TurboRegDialog.AFFINE: {
                targetJacobian = (targetPoint[1][0] - targetPoint[2][0])
                        * targetPoint[0][1]
                        + (targetPoint[2][0] - targetPoint[0][0])
                        * targetPoint[1][1]
                        + (targetPoint[0][0] - targetPoint[1][0])
                        * targetPoint[2][1];
                inverseMarquardtLevenbergOptimization(
                        iterationPower * maxIterations - 1);
                break;
            }
            case TurboRegDialog.BILINEAR: {
                marquardtLevenbergOptimization(
                        iterationPower * maxIterations - 1);
                break;
            }
        }
        scaleUpLandmarks();
        sourcePh.setPoints(sourcePoint);
        iterationCost *= ITERATION_PROGRESSION;
    }
    iterationPower /= ITERATION_PROGRESSION;
    if (transformation == TurboRegDialog.BILINEAR) {
        inNx = sourceImg.getWidth();
        inNy = sourceImg.getHeight();
        inImg = sourceImg.getCoefficient();
        if (sourceMsk == null) {
            inMsk = null;
        }
        else {
            inMsk = sourceMsk.getMask();
        }
        outNx = targetImg.getWidth();
        outNy = targetImg.getHeight();
        outImg = targetImg.getImage();
        outMsk = targetMsk.getMask();
    }
    else {
        inNx = targetImg.getWidth();
        inNy = targetImg.getHeight();
        inImg = targetImg.getCoefficient();
        inMsk = targetMsk.getMask();
        outNx = sourceImg.getWidth();
        outNy = sourceImg.getHeight();
        outImg = sourceImg.getImage();
        xGradient = sourceImg.getXGradient();
        yGradient = sourceImg.getYGradient();
        if (sourceMsk == null) {
            outMsk = null;
        }
        else {
            outMsk = sourceMsk.getMask();
        }
    }
    twiceInNx = 2 * inNx;
    twiceInNy = 2 * inNy;
    if (accelerated) {
        TurboRegProgressBar.skipProgressBar(
                iterationCost * (maxIterations - 1));
    }
    else {
        switch (transformation) {
            case TurboRegDialog.RIGID_BODY: {
                inverseMarquardtLevenbergRigidBodyOptimization(
                        maxIterations - 1);
                break;
            }
            case TurboRegDialog.TRANSLATION:
            case TurboRegDialog.SCALED_ROTATION:
            case TurboRegDialog.AFFINE: {
                inverseMarquardtLevenbergOptimization(maxIterations - 1);
                break;
            }
            case TurboRegDialog.BILINEAR: {
                marquardtLevenbergOptimization(maxIterations - 1);
                break;
            }
        }
    }
    sourcePh.setPoints(sourcePoint);
    iterationPower = (int)Math.pow(
            (double)ITERATION_PROGRESSION, (double)pyramidDepth);
    TurboRegProgressBar.workloadDone(
            pyramidDepth * maxIterations * iterationPower / ITERATION_PROGRESSION
                    - (iterationPower - 1) / (ITERATION_PROGRESSION - 1));
} /* end doRegistration */

/*********************************************************************
 Save the current landmarks into a text file and return the path
    and name of the file. Rigid format.
    @see TurboRegDialog#loadLandmarks()
    ********************************************************************/
String TurboRegTransform::saveTransformation (
        String filename
) {
    inNx = sourceImg.getWidth();
    inNy = sourceImg.getHeight();
    outNx = targetImg.getWidth();
    outNy = targetImg.getHeight();
    String path = "";
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
            case TurboRegDialog.TRANSLATION: {
                fw.write("TRANSLATION\n");
                break;
            }
            case TurboRegDialog.RIGID_BODY: {
                fw.write("RIGID_BODY\n");
                break;
            }
            case TurboRegDialog.SCALED_ROTATION: {
                fw.write("SCALED_ROTATION\n");
                break;
            }
            case TurboRegDialog.AFFINE: {
                fw.write("AFFINE\n");
                break;
            }
            case TurboRegDialog.BILINEAR: {
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
        if (transformation == TurboRegDialog.RIGID_BODY) {
            for (int i = 0; (i < transformation); i++) {
                fw.write(sourcePoint[i][0] + "\t" + sourcePoint[i][1] + "\n");
            }
        }
        else {
            for (int i = 0; (i < (transformation / 2)); i++) {
                fw.write(sourcePoint[i][0] + "\t" + sourcePoint[i][1] + "\n");
            }
        }
        fw.write("\n");
        fw.write("Target landmarks\n");
        if (transformation == TurboRegDialog.RIGID_BODY) {
            for (int i = 0; (i < transformation); i++) {
                fw.write(targetPoint[i][0] + "\t" + targetPoint[i][1] + "\n");
            }
        }
        else {
            for (int i = 0; (i < (transformation / 2)); i++) {
                fw.write(targetPoint[i][0] + "\t" + targetPoint[i][1] + "\n");
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
} /* end saveTransformation */

/*....................................................................
private methods
....................................................................*/
/*------------------------------------------------------------------*/
TurboRegTransform::affineTransform (
        double[][] matrix
) {
    double yx;
    double yy;
    double x0;
    double y0;
    int xMsk;
    int yMsk;
    int k = 0;
    TurboRegProgressBar.addWorkload(outNy);
    yx = matrix[0][0];
    yy = matrix[1][0];
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
            x0 += matrix[0][1];
            y0 += matrix[1][1];
        }
        yx += matrix[0][2];
        yy += matrix[1][2];
        TurboRegProgressBar.stepProgressBar();
    }
    TurboRegProgressBar.workloadDone(outNy);
} /* affineTransform */

/*------------------------------------------------------------------*/
TurboRegTransform::affineTransform (
        double[][] matrix,
        float[] outMsk
) {
    double yx;
    double yy;
    double x0;
    double y0;
    int xMsk;
    int yMsk;
    int k = 0;
    TurboRegProgressBar.addWorkload(outNy);
    yx = matrix[0][0];
    yy = matrix[1][0];
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
            x0 += matrix[0][1];
            y0 += matrix[1][1];
        }
        yx += matrix[0][2];
        yy += matrix[1][2];
        TurboRegProgressBar.stepProgressBar();
    }
    TurboRegProgressBar.workloadDone(outNy);
} /* affineTransform */

/*------------------------------------------------------------------*/
TurboRegTransform::bilinearTransform (
        double[][] matrix
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
    TurboRegProgressBar.addWorkload(outNy);
    yx = matrix[0][0];
    yy = matrix[1][0];
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
            x0 += matrix[0][1] + yxy;
            y0 += matrix[1][1] + yyy;
        }
        yx += matrix[0][2];
        yy += matrix[1][2];
        yxy += matrix[0][3];
        yyy += matrix[1][3];
        TurboRegProgressBar.stepProgressBar();
    }
    TurboRegProgressBar.workloadDone(outNy);
} /* bilinearTransform */

/*------------------------------------------------------------------*/
TurboRegTransform::bilinearTransform (
        double[][] matrix,
        float[] outMsk
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
    TurboRegProgressBar.addWorkload(outNy);
    yx = matrix[0][0];
    yy = matrix[1][0];
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
            x0 += matrix[0][1] + yxy;
            y0 += matrix[1][1] + yyy;
        }
        yx += matrix[0][2];
        yy += matrix[1][2];
        yxy += matrix[0][3];
        yyy += matrix[1][3];
        TurboRegProgressBar.stepProgressBar();
    }
    TurboRegProgressBar.workloadDone(outNy);
} /* bilinearTransform */

/*------------------------------------------------------------------*/
TurboRegTransform::computeBilinearGradientConstants (
) {
    double u1 = targetPoint[0][0];
    double u2 = targetPoint[1][0];
    double u3 = targetPoint[2][0];
    double u4 = targetPoint[3][0];
    double v1 = targetPoint[0][1];
    double v2 = targetPoint[1][1];
    double v3 = targetPoint[2][1];
    double v4 = targetPoint[3][1];
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
        double[][] sourcePoint,
        double[][] matrix
) {
    double u1 = sourcePoint[0][0];
    double u2 = sourcePoint[1][0];
    double u3 = sourcePoint[2][0];
    double v1 = sourcePoint[0][1];
    double v2 = sourcePoint[1][1];
    double v3 = sourcePoint[2][1];
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
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    return(meanSquares / ((double)area * Math.abs(det / targetJacobian)));
} /* getAffineMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getAffineMeanSquares (
        double[][] sourcePoint,
        double[][] matrix,
        double[] gradient
) {
    double u1 = sourcePoint[0][0];
    double u2 = sourcePoint[1][0];
    double u3 = sourcePoint[2][0];
    double v1 = sourcePoint[0][1];
    double v2 = sourcePoint[1][1];
    double v3 = sourcePoint[2][1];
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
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    return(meanSquares / ((double)area * Math.abs(det / targetJacobian)));
} /* getAffineMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getAffineMeanSquares (
        double[][] sourcePoint,
        double[][] matrix,
        double[][] hessian,
        double[] gradient
) {
    double u1 = sourcePoint[0][0];
    double u2 = sourcePoint[1][0];
    double u3 = sourcePoint[2][0];
    double v1 = sourcePoint[0][1];
    double v2 = sourcePoint[1][1];
    double v3 = sourcePoint[2][1];
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
            hessian[i][j] = 0.0;
        }
    }
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dx0 * dx0;
                        hessian[0][1] += dx0 * dy0;
                        hessian[0][2] += dx0 * dx1;
                        hessian[0][3] += dx0 * dy1;
                        hessian[0][4] += dx0 * dx2;
                        hessian[0][5] += dx0 * dy2;
                        hessian[1][1] += dy0 * dy0;
                        hessian[1][2] += dy0 * dx1;
                        hessian[1][3] += dy0 * dy1;
                        hessian[1][4] += dy0 * dx2;
                        hessian[1][5] += dy0 * dy2;
                        hessian[2][2] += dx1 * dx1;
                        hessian[2][3] += dx1 * dy1;
                        hessian[2][4] += dx1 * dx2;
                        hessian[2][5] += dx1 * dy2;
                        hessian[3][3] += dy1 * dy1;
                        hessian[3][4] += dy1 * dx2;
                        hessian[3][5] += dy1 * dy2;
                        hessian[4][4] += dx2 * dx2;
                        hessian[4][5] += dx2 * dy2;
                        hessian[5][5] += dy2 * dy2;
                    }
                }
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dx0 * dx0;
                        hessian[0][1] += dx0 * dy0;
                        hessian[0][2] += dx0 * dx1;
                        hessian[0][3] += dx0 * dy1;
                        hessian[0][4] += dx0 * dx2;
                        hessian[0][5] += dx0 * dy2;
                        hessian[1][1] += dy0 * dy0;
                        hessian[1][2] += dy0 * dx1;
                        hessian[1][3] += dy0 * dy1;
                        hessian[1][4] += dy0 * dx2;
                        hessian[1][5] += dy0 * dy2;
                        hessian[2][2] += dx1 * dx1;
                        hessian[2][3] += dx1 * dy1;
                        hessian[2][4] += dx1 * dx2;
                        hessian[2][5] += dx1 * dy2;
                        hessian[3][3] += dy1 * dy1;
                        hessian[3][4] += dy1 * dx2;
                        hessian[3][5] += dy1 * dy2;
                        hessian[4][4] += dx2 * dx2;
                        hessian[4][5] += dx2 * dy2;
                        hessian[5][5] += dy2 * dy2;
                    }
                }
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian[i][j] = hessian[j][i];
        }
    }
    return(meanSquares / ((double)area * Math.abs(det / targetJacobian)));
} /* getAffineMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getBilinearMeanSquares (
        double[][] matrix
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
    if (inMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1] + yxy;
                y0 += matrix[1][1] + yyy;
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
            yxy += matrix[0][3];
            yyy += matrix[1][3];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1] + yxy;
                y0 += matrix[1][1] + yyy;
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
            yxy += matrix[0][3];
            yyy += matrix[1][3];
        }
    }
    return(meanSquares / (double)area);
} /* getBilinearMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getBilinearMeanSquares (
        double[][] matrix,
        double[][] hessian,
        double[] gradient
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
            hessian[i][j] = 0.0;
        }
    }
    if (inMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dx0 * dx0;
                        hessian[0][1] += dx0 * dy0;
                        hessian[0][2] += dx0 * dx1;
                        hessian[0][3] += dx0 * dy1;
                        hessian[0][4] += dx0 * dx2;
                        hessian[0][5] += dx0 * dy2;
                        hessian[0][6] += dx0 * dx3;
                        hessian[0][7] += dx0 * dy3;
                        hessian[1][1] += dy0 * dy0;
                        hessian[1][2] += dy0 * dx1;
                        hessian[1][3] += dy0 * dy1;
                        hessian[1][4] += dy0 * dx2;
                        hessian[1][5] += dy0 * dy2;
                        hessian[1][6] += dy0 * dx3;
                        hessian[1][7] += dy0 * dy3;
                        hessian[2][2] += dx1 * dx1;
                        hessian[2][3] += dx1 * dy1;
                        hessian[2][4] += dx1 * dx2;
                        hessian[2][5] += dx1 * dy2;
                        hessian[2][6] += dx1 * dx3;
                        hessian[2][7] += dx1 * dy3;
                        hessian[3][3] += dy1 * dy1;
                        hessian[3][4] += dy1 * dx2;
                        hessian[3][5] += dy1 * dy2;
                        hessian[3][6] += dy1 * dx3;
                        hessian[3][7] += dy1 * dy3;
                        hessian[4][4] += dx2 * dx2;
                        hessian[4][5] += dx2 * dy2;
                        hessian[4][6] += dx2 * dx3;
                        hessian[4][7] += dx2 * dy3;
                        hessian[5][5] += dy2 * dy2;
                        hessian[5][6] += dy2 * dx3;
                        hessian[5][7] += dy2 * dy3;
                        hessian[6][6] += dx3 * dx3;
                        hessian[6][7] += dx3 * dy3;
                        hessian[7][7] += dy3 * dy3;
                    }
                }
                x0 += matrix[0][1] + yxy;
                y0 += matrix[1][1] + yyy;
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
            yxy += matrix[0][3];
            yyy += matrix[1][3];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dx0 * dx0;
                        hessian[0][1] += dx0 * dy0;
                        hessian[0][2] += dx0 * dx1;
                        hessian[0][3] += dx0 * dy1;
                        hessian[0][4] += dx0 * dx2;
                        hessian[0][5] += dx0 * dy2;
                        hessian[0][6] += dx0 * dx3;
                        hessian[0][7] += dx0 * dy3;
                        hessian[1][1] += dy0 * dy0;
                        hessian[1][2] += dy0 * dx1;
                        hessian[1][3] += dy0 * dy1;
                        hessian[1][4] += dy0 * dx2;
                        hessian[1][5] += dy0 * dy2;
                        hessian[1][6] += dy0 * dx3;
                        hessian[1][7] += dy0 * dy3;
                        hessian[2][2] += dx1 * dx1;
                        hessian[2][3] += dx1 * dy1;
                        hessian[2][4] += dx1 * dx2;
                        hessian[2][5] += dx1 * dy2;
                        hessian[2][6] += dx1 * dx3;
                        hessian[2][7] += dx1 * dy3;
                        hessian[3][3] += dy1 * dy1;
                        hessian[3][4] += dy1 * dx2;
                        hessian[3][5] += dy1 * dy2;
                        hessian[3][6] += dy1 * dx3;
                        hessian[3][7] += dy1 * dy3;
                        hessian[4][4] += dx2 * dx2;
                        hessian[4][5] += dx2 * dy2;
                        hessian[4][6] += dx2 * dx3;
                        hessian[4][7] += dx2 * dy3;
                        hessian[5][5] += dy2 * dy2;
                        hessian[5][6] += dy2 * dx3;
                        hessian[5][7] += dy2 * dy3;
                        hessian[6][6] += dx3 * dx3;
                        hessian[6][7] += dx3 * dy3;
                        hessian[7][7] += dy3 * dy3;
                    }
                }
                x0 += matrix[0][1] + yxy;
                y0 += matrix[1][1] + yyy;
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
            yxy += matrix[0][3];
            yyy += matrix[1][3];
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian[i][j] = hessian[j][i];
        }
    }
    return(meanSquares / (double)area);
} /* getBilinearMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getRigidBodyMeanSquares (
        double[][] matrix
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
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    return(meanSquares / (double)area);
} /* getRigidBodyMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getRigidBodyMeanSquares (
        double[][] matrix,
        double[] gradient

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
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    return(meanSquares / (double)area);
} /* getRigidBodyMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getRigidBodyMeanSquares (
        double[][] matrix,
        double[][] hessian,
        double[] gradient
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
            hessian[i][j] = 0.0;
        }
    }
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dTheta * dTheta;
                        hessian[0][1] += dTheta * xGradient[k];
                        hessian[0][2] += dTheta * yGradient[k];
                        hessian[1][1] += xGradient[k] * xGradient[k];
                        hessian[1][2] += xGradient[k] * yGradient[k];
                        hessian[2][2] += yGradient[k] * yGradient[k];
                    }
                }
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dTheta * dTheta;
                        hessian[0][1] += dTheta * xGradient[k];
                        hessian[0][2] += dTheta * yGradient[k];
                        hessian[1][1] += xGradient[k] * xGradient[k];
                        hessian[1][2] += xGradient[k] * yGradient[k];
                        hessian[2][2] += yGradient[k] * yGradient[k];
                    }
                }
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian[i][j] = hessian[j][i];
        }
    }
    return(meanSquares / (double)area);
} /* getRigidBodyMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getScaledRotationMeanSquares (
        double[][] sourcePoint,
        double[][] matrix
) {
    double u1 = sourcePoint[0][0];
    double u2 = sourcePoint[1][0];
    double v1 = sourcePoint[0][1];
    double v2 = sourcePoint[1][1];
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
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    return(meanSquares / ((double)area * uv2 / targetJacobian));
} /* getScaledRotationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getScaledRotationMeanSquares (
        double[][] sourcePoint,
        double[][] matrix,
        double[] gradient
) {
    double u1 = sourcePoint[0][0];
    double u2 = sourcePoint[1][0];
    double v1 = sourcePoint[0][1];
    double v2 = sourcePoint[1][1];
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
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    return(meanSquares / ((double)area * uv2 / targetJacobian));
} /* getScaledRotationMeanSquares */

/*------------------------------------------------------------------*/
double TurboRegTransform::getScaledRotationMeanSquares (
        double[][] sourcePoint,
        double[][] matrix,
        double[][] hessian,
        double[] gradient
) {
    double u1 = sourcePoint[0][0];
    double u2 = sourcePoint[1][0];
    double v1 = sourcePoint[0][1];
    double v2 = sourcePoint[1][1];
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
            hessian[i][j] = 0.0;
        }
    }
    if (outMsk == null) {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dx0 * dx0;
                        hessian[0][1] += dx0 * dy0;
                        hessian[0][2] += dx0 * dx1;
                        hessian[0][3] += dx0 * dy1;
                        hessian[1][1] += dy0 * dy0;
                        hessian[1][2] += dy0 * dx1;
                        hessian[1][3] += dy0 * dy1;
                        hessian[2][2] += dx1 * dx1;
                        hessian[2][3] += dx1 * dy1;
                        hessian[3][3] += dy1 * dy1;
                    }
                }
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    else {
        yx = matrix[0][0];
        yy = matrix[1][0];
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
                        hessian[0][0] += dx0 * dx0;
                        hessian[0][1] += dx0 * dy0;
                        hessian[0][2] += dx0 * dx1;
                        hessian[0][3] += dx0 * dy1;
                        hessian[1][1] += dy0 * dy0;
                        hessian[1][2] += dy0 * dx1;
                        hessian[1][3] += dy0 * dy1;
                        hessian[2][2] += dx1 * dx1;
                        hessian[2][3] += dx1 * dy1;
                        hessian[3][3] += dy1 * dy1;
                    }
                }
                x0 += matrix[0][1];
                y0 += matrix[1][1];
            }
            yx += matrix[0][2];
            yy += matrix[1][2];
        }
    }
    for (int i = 1; (i < transformation); i++) {
        for (int j = 0; (j < i); j++) {
            hessian[i][j] = hessian[j][i];
        }
    }
    return(meanSquares / ((double)area * uv2 / targetJacobian));
} /* getScaledRotationMeanSquares */

/*------------------------------------------------------------------*/
private double[][] getTransformationMatrix (
        double[][] fromCoord,
        double[][] toCoord
) {
    double[][] matrix = null;
    double[][] a = null;
    double[] v = null;
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            matrix = new double[2][1];
            matrix[0][0] = toCoord[0][0] - fromCoord[0][0];
            matrix[1][0] = toCoord[0][1] - fromCoord[0][1];
            break;
        }
        case TurboRegDialog.RIGID_BODY: {
            double angle = Math.atan2(fromCoord[2][0] - fromCoord[1][0],
                    fromCoord[2][1] - fromCoord[1][1])
                    - Math.atan2(toCoord[2][0] - toCoord[1][0],
                    toCoord[2][1] - toCoord[1][1]);
            double c = Math.cos(angle);
            double s = Math.sin(angle);
            matrix = new double[2][3];
            matrix[0][0] = toCoord[0][0]
                    - c * fromCoord[0][0] + s * fromCoord[0][1];
            matrix[0][1] = c;
            matrix[0][2] = -s;
            matrix[1][0] = toCoord[0][1]
                    - s * fromCoord[0][0] - c * fromCoord[0][1];
            matrix[1][1] = s;
            matrix[1][2] = c;
            break;
        }
        case TurboRegDialog.SCALED_ROTATION: {
            matrix = new double[2][3];
            a = new double[3][3];
            v = new double[3];
            a[0][0] = 1.0;
            a[0][1] = fromCoord[0][0];
            a[0][2] = fromCoord[0][1];
            a[1][0] = 1.0;
            a[1][1] = fromCoord[1][0];
            a[1][2] = fromCoord[1][1];
            a[2][0] = 1.0;
            a[2][1] = fromCoord[0][1] - fromCoord[1][1] + fromCoord[1][0];
            a[2][2] = fromCoord[1][0] + fromCoord[1][1] - fromCoord[0][0];
            invertGauss(a);
            v[0] = toCoord[0][0];
            v[1] = toCoord[1][0];
            v[2] = toCoord[0][1] - toCoord[1][1] + toCoord[1][0];
            for (int i = 0; (i < 3); i++) {
                matrix[0][i] = 0.0;
                for (int j = 0; (j < 3); j++) {
                    matrix[0][i] += a[i][j] * v[j];
                }
            }
            v[0] = toCoord[0][1];
            v[1] = toCoord[1][1];
            v[2] = toCoord[1][0] + toCoord[1][1] - toCoord[0][0];
            for (int i = 0; (i < 3); i++) {
                matrix[1][i] = 0.0;
                for (int j = 0; (j < 3); j++) {
                    matrix[1][i] += a[i][j] * v[j];
                }
            }
            break;
        }
        case TurboRegDialog.AFFINE: {
            matrix = new double[2][3];
            a = new double[3][3];
            v = new double[3];
            a[0][0] = 1.0;
            a[0][1] = fromCoord[0][0];
            a[0][2] = fromCoord[0][1];
            a[1][0] = 1.0;
            a[1][1] = fromCoord[1][0];
            a[1][2] = fromCoord[1][1];
            a[2][0] = 1.0;
            a[2][1] = fromCoord[2][0];
            a[2][2] = fromCoord[2][1];
            invertGauss(a);
            v[0] = toCoord[0][0];
            v[1] = toCoord[1][0];
            v[2] = toCoord[2][0];
            for (int i = 0; (i < 3); i++) {
                matrix[0][i] = 0.0;
                for (int j = 0; (j < 3); j++) {
                    matrix[0][i] += a[i][j] * v[j];
                }
            }
            v[0] = toCoord[0][1];
            v[1] = toCoord[1][1];
            v[2] = toCoord[2][1];
            for (int i = 0; (i < 3); i++) {
                matrix[1][i] = 0.0;
                for (int j = 0; (j < 3); j++) {
                    matrix[1][i] += a[i][j] * v[j];
                }
            }
            break;
        }
        case TurboRegDialog.BILINEAR: {
            matrix = new double[2][4];
            a = new double[4][4];
            v = new double[4];
            a[0][0] = 1.0;
            a[0][1] = fromCoord[0][0];
            a[0][2] = fromCoord[0][1];
            a[0][3] = fromCoord[0][0] * fromCoord[0][1];
            a[1][0] = 1.0;
            a[1][1] = fromCoord[1][0];
            a[1][2] = fromCoord[1][1];
            a[1][3] = fromCoord[1][0] * fromCoord[1][1];
            a[2][0] = 1.0;
            a[2][1] = fromCoord[2][0];
            a[2][2] = fromCoord[2][1];
            a[2][3] = fromCoord[2][0] * fromCoord[2][1];
            a[3][0] = 1.0;
            a[3][1] = fromCoord[3][0];
            a[3][2] = fromCoord[3][1];
            a[3][3] = fromCoord[3][0] * fromCoord[3][1];
            invertGauss(a);
            v[0] = toCoord[0][0];
            v[1] = toCoord[1][0];
            v[2] = toCoord[2][0];
            v[3] = toCoord[3][0];
            for (int i = 0; (i < 4); i++) {
                matrix[0][i] = 0.0;
                for (int j = 0; (j < 4); j++) {
                    matrix[0][i] += a[i][j] * v[j];
                }
            }
            v[0] = toCoord[0][1];
            v[1] = toCoord[1][1];
            v[2] = toCoord[2][1];
            v[3] = toCoord[3][1];
            for (int i = 0; (i < 4); i++) {
                matrix[1][i] = 0.0;
                for (int j = 0; (j < 4); j++) {
                    matrix[1][i] += a[i][j] * v[j];
                }
            }
            break;
        }
    }
    return(matrix);
} /* end getTransformationMatrix */

/*------------------------------------------------------------------*/
double TurboRegTransform::getTranslationMeanSquares (
        double[][] matrix
) {
    double dx = matrix[0][0];
    double dy = matrix[1][0];
    double dx0 = dx;
    double difference;
    double meanSquares = 0.0;
    long area = 0L;
    int xMsk;
    int yMsk;
    int k = 0;
    x = dx - Math.floor(dx);
    y = dy - Math.floor(dy);
    xWeights();
    yWeights();
    if (outMsk == null) {
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
        double[][] matrix,
        double[] gradient
) {
    double dx = matrix[0][0];
    double dy = matrix[1][0];
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
    x = dx - Math.floor(dx);
    y = dy - Math.floor(dy);
    xWeights();
    yWeights();
    if (outMsk == null) {
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
        double[][] matrix,
        double[][] hessian,
        double[] gradient
) {
    double dx = matrix[0][0];
    double dy = matrix[1][0];
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
            hessian[i][j] = 0.0;
        }
    }
    x = dx - Math.floor(dx);
    y = dy - Math.floor(dy);
    xWeights();
    yWeights();
    if (outMsk == null) {
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
                            hessian[0][0] += xGradient[k] * xGradient[k];
                            hessian[0][1] += xGradient[k] * yGradient[k];
                            hessian[1][1] += yGradient[k] * yGradient[k];
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
                            hessian[0][0] += xGradient[k] * xGradient[k];
                            hessian[0][1] += xGradient[k] * yGradient[k];
                            hessian[1][1] += yGradient[k] * yGradient[k];
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
            hessian[i][j] = hessian[j][i];
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
TurboRegTransform::inverseMarquardtLevenbergOptimization (
        int workload
) {
    double[][] attempt = new double[transformation / 2][2];
    double[][] hessian = new double[transformation][transformation];
    double[][] pseudoHessian = new double[transformation][transformation];
    double[] gradient = new double[transformation];
    double[][] matrix = getTransformationMatrix(sourcePoint, targetPoint);
    double[] update = new double[transformation];
    double bestMeanSquares = 0.0;
    double meanSquares = 0.0;
    double lambda = FIRST_LAMBDA;
    double displacement;
    int iteration = 0;
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            bestMeanSquares = getTranslationMeanSquares(
                    matrix, hessian, gradient);
            break;
        }
        case TurboRegDialog.SCALED_ROTATION: {
            bestMeanSquares = getScaledRotationMeanSquares(
                    sourcePoint, matrix, hessian, gradient);
            break;
        }
        case TurboRegDialog.AFFINE: {
            bestMeanSquares = getAffineMeanSquares(
                    sourcePoint, matrix, hessian, gradient);
            break;
        }
    }
    iteration++;
    do {
        for (int k = 0; (k < transformation); k++) {
            pseudoHessian[k][k] = (1.0 + lambda) * hessian[k][k];
        }
        invertGauss(pseudoHessian);
        update = matrixMultiply(pseudoHessian, gradient);
        displacement = 0.0;
        for (int k = 0; (k < (transformation / 2)); k++) {
            attempt[k][0] = sourcePoint[k][0] - update[2 * k];
            attempt[k][1] = sourcePoint[k][1] - update[2 * k + 1];
            displacement += Math.sqrt(update[2 * k] * update[2 * k]
                    + update[2 * k + 1] * update[2 * k + 1]);
        }
        displacement /= 0.5 * (double)transformation;
        matrix = getTransformationMatrix(attempt, targetPoint);
        switch (transformation) {
            case TurboRegDialog.TRANSLATION: {
                if (accelerated) {
                    meanSquares = getTranslationMeanSquares(
                            matrix, gradient);
                }
                else {
                    meanSquares = getTranslationMeanSquares(
                            matrix, hessian, gradient);
                }
                break;
            }
            case TurboRegDialog.SCALED_ROTATION: {
                if (accelerated) {
                    meanSquares = getScaledRotationMeanSquares(
                            attempt, matrix, gradient);
                }
                else {
                    meanSquares = getScaledRotationMeanSquares(
                            attempt, matrix, hessian, gradient);
                }
                break;
            }
            case TurboRegDialog.AFFINE: {
                if (accelerated) {
                    meanSquares = getAffineMeanSquares(
                            attempt, matrix, gradient);
                }
                else {
                    meanSquares = getAffineMeanSquares(
                            attempt, matrix, hessian, gradient);
                }
                break;
            }
        }
        iteration++;
        if (meanSquares < bestMeanSquares) {
            bestMeanSquares = meanSquares;
            for (int k = 0; (k < (transformation / 2)); k++) {
                sourcePoint[k][0] = attempt[k][0];
                sourcePoint[k][1] = attempt[k][1];
            }
            lambda /= LAMBDA_MAGSTEP;
        }
        else {
            lambda *= LAMBDA_MAGSTEP;
        }
        TurboRegProgressBar.skipProgressBar(iterationCost);
        workload--;
    } while ((iteration < (maxIterations * iterationPower - 1))
            && (pixelPrecision <= displacement));
    invertGauss(hessian);
    update = matrixMultiply(hessian, gradient);
    for (int k = 0; (k < (transformation / 2)); k++) {
        attempt[k][0] = sourcePoint[k][0] - update[2 * k];
        attempt[k][1] = sourcePoint[k][1] - update[2 * k + 1];
    }
    matrix = getTransformationMatrix(attempt, targetPoint);
    switch (transformation) {
        case TurboRegDialog.TRANSLATION: {
            meanSquares = getTranslationMeanSquares(matrix);
            break;
        }
        case TurboRegDialog.SCALED_ROTATION: {
            meanSquares = getScaledRotationMeanSquares(attempt, matrix);
            break;
        }
        case TurboRegDialog.AFFINE: {
            meanSquares = getAffineMeanSquares(attempt, matrix);
            break;
        }
    }
    iteration++;
    if (meanSquares < bestMeanSquares) {
        for (int k = 0; (k < (transformation / 2)); k++) {
            sourcePoint[k][0] = attempt[k][0];
            sourcePoint[k][1] = attempt[k][1];
        }
    }
    TurboRegProgressBar.skipProgressBar(workload * iterationCost);
} /* end inverseMarquardtLevenbergOptimization */

/*------------------------------------------------------------------*/
TurboRegTransform::inverseMarquardtLevenbergRigidBodyOptimization (
        int workload
) {
    double[][] attempt = new double[2][3];
    double[][] hessian = new double[transformation][transformation];
    double[][] pseudoHessian = new double[transformation][transformation];
    double[] gradient = new double[transformation];
    double[][] matrix = getTransformationMatrix(targetPoint, sourcePoint);
    double[] update = new double[transformation];
    double bestMeanSquares = 0.0;
    double meanSquares = 0.0;
    double lambda = FIRST_LAMBDA;
    double angle;
    double c;
    double s;
    double displacement;
    int iteration = 0;
    for (int k = 0; (k < transformation); k++) {
        sourcePoint[k][0] = matrix[0][0] + targetPoint[k][0] * matrix[0][1]
                + targetPoint[k][1] * matrix[0][2];
        sourcePoint[k][1] = matrix[1][0] + targetPoint[k][0] * matrix[1][1]
                + targetPoint[k][1] * matrix[1][2];
    }
    matrix = getTransformationMatrix(sourcePoint, targetPoint);
    bestMeanSquares = getRigidBodyMeanSquares(matrix, hessian, gradient);
    iteration++;
    do {
        for (int k = 0; (k < transformation); k++) {
            pseudoHessian[k][k] = (1.0 + lambda) * hessian[k][k];
        }
        invertGauss(pseudoHessian);
        update = matrixMultiply(pseudoHessian, gradient);
        angle = Math.atan2(matrix[0][2], matrix[0][1]) - update[0];
        attempt[0][1] = Math.cos(angle);
        attempt[0][2] = Math.sin(angle);
        attempt[1][1] = -attempt[0][2];
        attempt[1][2] = attempt[0][1];
        c = Math.cos(update[0]);
        s = Math.sin(update[0]);
        attempt[0][0] = (matrix[0][0] + update[1]) * c
                - (matrix[1][0] + update[2]) * s;
        attempt[1][0] = (matrix[0][0] + update[1]) * s
                + (matrix[1][0] + update[2]) * c;
        displacement = Math.sqrt(update[1] * update[1] + update[2] * update[2])
                + 0.25 * Math.sqrt((double)(inNx * inNx) + (double)(inNy * inNy))
                * Math.abs(update[0]);
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
                    matrix[i][j] = attempt[i][j];
                }
            }
            lambda /= LAMBDA_MAGSTEP;
        }
        else {
            lambda *= LAMBDA_MAGSTEP;
        }
        TurboRegProgressBar.skipProgressBar(iterationCost);
        workload--;
    } while ((iteration < (maxIterations * iterationPower - 1))
            && (pixelPrecision <= displacement));
    invertGauss(hessian);
    update = matrixMultiply(hessian, gradient);
    angle = Math.atan2(matrix[0][2], matrix[0][1]) - update[0];
    attempt[0][1] = Math.cos(angle);
    attempt[0][2] = Math.sin(angle);
    attempt[1][1] = -attempt[0][2];
    attempt[1][2] = attempt[0][1];
    c = Math.cos(update[0]);
    s = Math.sin(update[0]);
    attempt[0][0] = (matrix[0][0] + update[1]) * c
            - (matrix[1][0] + update[2]) * s;
    attempt[1][0] = (matrix[0][0] + update[1]) * s
            + (matrix[1][0] + update[2]) * c;
    meanSquares = getRigidBodyMeanSquares(attempt);
    iteration++;
    if (meanSquares < bestMeanSquares) {
        for (int i = 0; (i < 2); i++) {
            for (int j = 0; (j < 3); j++) {
                matrix[i][j] = attempt[i][j];
            }
        }
    }
    for (int k = 0; (k < transformation); k++) {
        sourcePoint[k][0] = (targetPoint[k][0] - matrix[0][0]) * matrix[0][1]
                + (targetPoint[k][1] - matrix[1][0]) * matrix[1][1];
        sourcePoint[k][1] = (targetPoint[k][0] - matrix[0][0]) * matrix[0][2]
                + (targetPoint[k][1] - matrix[1][0]) * matrix[1][2];
    }
    TurboRegProgressBar.skipProgressBar(workload * iterationCost);
} /* end inverseMarquardtLevenbergRigidBodyOptimization */

/*------------------------------------------------------------------*/
TurboRegTransform::invertGauss (
        double[][] matrix
) {
    int n = matrix.length;
    double[][] inverse = new double[n][n];
    for (int i = 0; (i < n); i++) {
        double max = matrix[i][0];
        double absMax = Math.abs(max);
        for (int j = 0; (j < n); j++) {
            inverse[i][j] = 0.0;
            if (absMax < Math.abs(matrix[i][j])) {
                max = matrix[i][j];
                absMax = Math.abs(max);
            }
        }
        inverse[i][i] = 1.0 / max;
        for (int j = 0; (j < n); j++) {
            matrix[i][j] /= max;
        }
    }
    for (int j = 0; (j < n); j++) {
        double max = matrix[j][j];
        double absMax = Math.abs(max);
        int k = j;
        for (int i = j + 1; (i < n); i++) {
            if (absMax < Math.abs(matrix[i][j])) {
                max = matrix[i][j];
                absMax = Math.abs(max);
                k = i;
            }
        }
        if (k != j) {
            double[] partialLine = new double[n - j];
            double[] fullLine = new double[n];
            System.arraycopy(matrix[j], j, partialLine, 0, n - j);
            System.arraycopy(matrix[k], j, matrix[j], j, n - j);
            System.arraycopy(partialLine, 0, matrix[k], j, n - j);
            System.arraycopy(inverse[j], 0, fullLine, 0, n);
            System.arraycopy(inverse[k], 0, inverse[j], 0, n);
            System.arraycopy(fullLine, 0, inverse[k], 0, n);
        }
        for (k = 0; (k <= j); k++) {
            inverse[j][k] /= max;
        }
        for (k = j + 1; (k < n); k++) {
            matrix[j][k] /= max;
            inverse[j][k] /= max;
        }
        for (int i = j + 1; (i < n); i++) {
            for (k = 0; (k <= j); k++) {
                inverse[i][k] -= matrix[i][j] * inverse[j][k];
            }
            for (k = j + 1; (k < n); k++) {
                matrix[i][k] -= matrix[i][j] * matrix[j][k];
                inverse[i][k] -= matrix[i][j] * inverse[j][k];
            }
        }
    }
    for (int j = n - 1; (1 <= j); j--) {
        for (int i = j - 1; (0 <= i); i--) {
            for (int k = 0; (k <= j); k++) {
                inverse[i][k] -= matrix[i][j] * inverse[j][k];
            }
            for (int k = j + 1; (k < n); k++) {
                matrix[i][k] -= matrix[i][j] * matrix[j][k];
                inverse[i][k] -= matrix[i][j] * inverse[j][k];
            }
        }
    }
    for (int i = 0; (i < n); i++) {
        System.arraycopy(inverse[i], 0, matrix[i], 0, n);
    }
} /* end invertGauss */

/*------------------------------------------------------------------*/
TurboRegTransform::marquardtLevenbergOptimization (
        int workload
) {
    double[][] attempt = new double[transformation / 2][2];
    double[][] hessian = new double[transformation][transformation];
    double[][] pseudoHessian = new double[transformation][transformation];
    double[] gradient = new double[transformation];
    double[][] matrix = getTransformationMatrix(targetPoint, sourcePoint);
    double[] update = new double[transformation];
    double bestMeanSquares = 0.0;
    double meanSquares = 0.0;
    double lambda = FIRST_LAMBDA;
    double displacement;
    int iteration = 0;
    bestMeanSquares = getBilinearMeanSquares(matrix, hessian, gradient);
    iteration++;
    do {
        for (int k = 0; (k < transformation); k++) {
            pseudoHessian[k][k] = (1.0 + lambda) * hessian[k][k];
        }
        invertGauss(pseudoHessian);
        update = matrixMultiply(pseudoHessian, gradient);
        displacement = 0.0;
        for (int k = 0; (k < (transformation / 2)); k++) {
            attempt[k][0] = sourcePoint[k][0] - update[2 * k];
            attempt[k][1] = sourcePoint[k][1] - update[2 * k + 1];
            displacement += Math.sqrt(update[2 * k] * update[2 * k]
                    + update[2 * k + 1] * update[2 * k + 1]);
        }
        displacement /= 0.5 * (double)transformation;
        matrix = getTransformationMatrix(targetPoint, attempt);
        meanSquares = getBilinearMeanSquares(matrix, hessian, gradient);
        iteration++;
        if (meanSquares < bestMeanSquares) {
            bestMeanSquares = meanSquares;
            for (int k = 0; (k < (transformation / 2)); k++) {
                sourcePoint[k][0] = attempt[k][0];
                sourcePoint[k][1] = attempt[k][1];
            }
            lambda /= LAMBDA_MAGSTEP;
        }
        else {
            lambda *= LAMBDA_MAGSTEP;
        }
        TurboRegProgressBar.skipProgressBar(iterationCost);
        workload--;
    } while ((iteration < (maxIterations * iterationPower - 1))
            && (pixelPrecision <= displacement));
    invertGauss(hessian);
    update = matrixMultiply(hessian, gradient);
    for (int k = 0; (k < (transformation / 2)); k++) {
        attempt[k][0] = sourcePoint[k][0] - update[2 * k];
        attempt[k][1] = sourcePoint[k][1] - update[2 * k + 1];
    }
    matrix = getTransformationMatrix(targetPoint, attempt);
    meanSquares = getBilinearMeanSquares(matrix);
    iteration++;
    if (meanSquares < bestMeanSquares) {
        for (int k = 0; (k < (transformation / 2)); k++) {
            sourcePoint[k][0] = attempt[k][0];
            sourcePoint[k][1] = attempt[k][1];
        }
    }
    TurboRegProgressBar.skipProgressBar(workload * iterationCost);
} /* end marquardtLevenbergOptimization */

/*------------------------------------------------------------------*/
private double[] matrixMultiply (
        double[][] matrix,
        double[] vector
) {
    double[] result = new double[matrix.length];
    for (int i = 0; (i < matrix.length); i++) {
        result[i] = 0.0;
        for (int j = 0; (j < vector.length); j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return(result);
} /* end matrixMultiply */

/*------------------------------------------------------------------*/
TurboRegTransform::scaleBottomDownLandmarks (
) {
    for (int depth = 1; (depth < pyramidDepth); depth++) {
        if (transformation == TurboRegDialog.RIGID_BODY) {
            for (int n = 0; (n < transformation); n++) {
                sourcePoint[n][0] *= 0.5;
                sourcePoint[n][1] *= 0.5;
                targetPoint[n][0] *= 0.5;
                targetPoint[n][1] *= 0.5;
            }
        }
        else {
            for (int n = 0; (n < (transformation / 2)); n++) {
                sourcePoint[n][0] *= 0.5;
                sourcePoint[n][1] *= 0.5;
                targetPoint[n][0] *= 0.5;
                targetPoint[n][1] *= 0.5;
            }
        }
    }
} /* end scaleBottomDownLandmarks */

/*------------------------------------------------------------------*/
TurboRegTransform::scaleUpLandmarks (
) {
    if (transformation == TurboRegDialog.RIGID_BODY) {
        for (int n = 0; (n < transformation); n++) {
            sourcePoint[n][0] *= 2.0;
            sourcePoint[n][1] *= 2.0;
            targetPoint[n][0] *= 2.0;
            targetPoint[n][1] *= 2.0;
        }
    }
    else {
        for (int n = 0; (n < (transformation / 2)); n++) {
            sourcePoint[n][0] *= 2.0;
            sourcePoint[n][1] *= 2.0;
            targetPoint[n][0] *= 2.0;
            targetPoint[n][1] *= 2.0;
        }
    }
} /* end scaleUpLandmarks */

/*------------------------------------------------------------------*/
TurboRegTransform::translationTransform (
        double[][] matrix
) {
    double dx = matrix[0][0];
    double dy = matrix[1][0];
    double dx0 = dx;
    int xMsk;
    int yMsk;
    x = dx - Math.floor(dx);
    y = dy - Math.floor(dy);
    if (!accelerated) {
        xWeights();
        yWeights();
    }
    int k = 0;
    TurboRegProgressBar.addWorkload(outNy);
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
        TurboRegProgressBar.stepProgressBar();
    }
    TurboRegProgressBar.workloadDone(outNy);
} /* translationTransform */

/*------------------------------------------------------------------*/
TurboRegTransform::translationTransform (
        double[][] matrix,
        float[] outMsk
) {
    double dx = matrix[0][0];
    double dy = matrix[1][0];
    double dx0 = dx;
    int xMsk;
    int yMsk;
    x = dx - Math.floor(dx);
    y = dy - Math.floor(dy);
    if (!accelerated) {
        xWeights();
        yWeights();
    }
    int k = 0;
    TurboRegProgressBar.addWorkload(outNy);
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
        TurboRegProgressBar.stepProgressBar();
    }
    TurboRegProgressBar.workloadDone(outNy);
} /* translationTransform */

/*------------------------------------------------------------------*/
TurboRegTransform::xDxWeights (
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
TurboRegTransform::xIndexes (
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
TurboRegTransform::xWeights (
) {
    s = 1.0 - x;
    xWeight[3] = s * s * s / 6.0;
    s = x * x;
    xWeight[2] = 2.0 / 3.0 - 0.5 * s * (2.0 - x);
    xWeight[0] = s * x / 6.0;
    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
} /* xWeights */

/*------------------------------------------------------------------*/
TurboRegTransform::yDyWeights (
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
TurboRegTransform::yIndexes (
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
TurboRegTransform::yWeights (
) {
    t = 1.0 - y;
    yWeight[3] = t * t * t / 6.0;
    t = y * y;
    yWeight[2] = 2.0 / 3.0 - 0.5 * t * (2.0 - y);
    yWeight[0] = t * y / 6.0;
    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
} /* yWeights */

