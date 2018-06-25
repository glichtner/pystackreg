#include "TurboReg.h"
#include "TurboRegImage.h"

#include <cstring>
#include <vector>
#include <cmath>

#define COUNTOF(x) (sizeof(x) / sizeof(*x))

/*====================================================================
|	turboRegImage
\===================================================================*/

/*********************************************************************
 This class is responsible for the image preprocessing that takes
place concurrently with user-interface events. It contains methods
to compute B-spline coefficients and their pyramids, image pyramids,
gradients, and gradient pyramids.
********************************************************************/


/*********************************************************************
 Start the image precomputations. The computation of the B-spline
    coefficients of the full-size image is not interruptible; all other
    methods are.
    ********************************************************************/
void TurboRegImage::init (
) {
    coefficient = getBasicFromCardinal2D();

    switch (transformation) {
        case GENERIC_TRANSFORMATION: {
            break;
        }
        case TRANSLATION:
        case RIGID_BODY:
        case SCALED_ROTATION:
        case AFFINE: {
            if (isTarget) {
                buildCoefficientPyramid();
            }
            else {
                imageToXYGradient2D();
                buildImageAndGradientPyramid();
            }
            break;
        }
        case BILINEAR: {
            if (isTarget) {
                buildImagePyramid();
            }
            else {
                buildCoefficientPyramid();
            }
            break;
        }
    }
} /* end run */

/*....................................................................
constructors
....................................................................*/
/*********************************************************************
 Converts the pixel array of the incoming <code>ImagePlus</code>
    object into a local <code>float</code> array.
    @param imp <code>ImagePlus</code> object to preprocess.
    @param transformation Transformation code.
    @param isTarget Tags the current object as a target or source image.
    ********************************************************************/
TurboRegImage::TurboRegImage (
    double *img, int width, int height, Transformation transformation, bool isTarget
) {
    transformation = transformation;
    isTarget = isTarget;
    width = width;
    height = height;
    int k = 0;

    image = new double[width * height];
    for (int y = 0; (y < height); y++) {
        for (int x = 0; (x < width); x++, k++) {
            image[k] = (double)img[k];
        }
    }

    init();
    
    /*if (imp.getType() == ImagePlus.GRAY8) {
        image = new float[width * height];
        byte[] pixels = (byte[])imp.getProcessor().getPixels();
        for (int y = 0; (y < height); y++) {
            for (int x = 0; (x < width); x++, k++) {
                image[k] = (float)(pixels[k] & 0xFF);
            }
            
        }
    }
    else if (imp.getType() == ImagePlus.GRAY16) {
        image = new float[width * height];
        short[] pixels = (short[])imp.getProcessor().getPixels();
        for (int y = 0; (y < height); y++) {
            for (int x = 0; (x < width); x++, k++) {
                if (pixels[k] < (short)0) {
                    image[k] = (float)pixels[k] + 65536.0F;
                }
                else {
                    image[k] = (float)pixels[k];
                }
            }
            
        }
    }
    else if (imp.getType() == ImagePlus.GRAY32) {
        image = (double*)imp.getProcessor().getPixels();
    }*/
    
} /* end turboRegImage */


void TurboRegImage::antiSymmetricFirMirrorOffBounds1D(
        double* h,
        int hlen,
        double* c,
        int clen,
        double* s,
        int slen
) {
    if (2 <= hlen) {
        s[0] = h[1] * (c[1] - c[0]);
        for (int i = 1; (i < (slen - 1)); i++) {
            s[i] = h[1] * (c[i + 1] - c[i - 1]);
        }
        s[slen - 1] = h[1] * (c[clen - 1] - c[clen - 2]);
    }
    else {
        s[0] = 0.0;
    }
} /* end antiSymmetricFirMirrorOffBounds1D */

/*------------------------------------------------------------------*/
void TurboRegImage::basicToCardinal2D (
        double* basic,
        double* cardinal,
        int width,
        int height,
        int degree
) {
    double* hLine = new double[width];
    double* vLine = new double[height];
    double* hData = new double[width];
    double* vData = new double[height];
    double* h = NULL;

    switch (degree) {
        case 3: {
            h = new double[2];
            h[0] = 2.0 / 3.0;
            h[1] = 1.0 / 6.0;
            break;
        }
        case 7: {
            h = new double[4];
            h[0] = 151.0 / 315.0;
            h[1] = 397.0 / 1680.0;
            h[2] = 1.0 / 42.0;
            h[3] = 1.0 / 5040.0;
            break;
        }
        default: {
            h = new double[1];
            h[0] = 1.0;
        }
    }
    
    for (int y = 0; (y < height); y++) {
        extractRow(basic, y, hLine, width);
        symmetricFirMirrorOffBounds1D(h, COUNTOF(h), hLine, width, hData, width);
        putRow(cardinal, y, hData, width);
    }

    for (int x = 0; (x < width); x++) {
        extractColumn(cardinal, width, x, vLine, height);
        symmetricFirMirrorOffBounds1D(h, COUNTOF(h), vLine, height, vData, height);
        putColumn(cardinal, width, x, vData, height);

    }

} /* end basicToCardinal2D */

/*------------------------------------------------------------------*/
void TurboRegImage::buildCoefficientPyramid (
) {
    int fullWidth;
    int fullHeight;
    double* fullDual = new double[width * height];
    int halfWidth = width;
    int halfHeight = height;
    if (1 < pyramidDepth) {
        basicToCardinal2D(coefficient, fullDual, width, height, 7);
    }
    for (int depth = 1; ((depth < pyramidDepth)); depth++) {
        fullWidth = halfWidth;
        fullHeight = halfHeight;
        halfWidth /= 2;
        halfHeight /= 2;
        double* halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
        double* halfCoefficient = getBasicFromCardinal2D(
                halfDual, halfWidth, halfHeight, 7);
        pyramid.push(halfCoefficient);
        pyramid.push(new Integer(halfHeight));
        pyramid.push(new Integer(halfWidth));
        fullDual = halfDual;
    }
} /* end buildCoefficientPyramid */

/*------------------------------------------------------------------*/
void TurboRegImage::buildImageAndGradientPyramid (
) {
    int fullWidth;
    int fullHeight;
    double* fullDual = new double[width * height];
    int halfWidth = width;
    int halfHeight = height;
    if (1 < pyramidDepth) {
        cardinalToDual2D(image, fullDual, width, height, 3);
    }
    for (int depth = 1; ((depth < pyramidDepth));
            depth++) {
        fullWidth = halfWidth;
        fullHeight = halfHeight;
        halfWidth /= 2;
        halfHeight /= 2;
        double* halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
        double* halfImage = getBasicFromCardinal2D(
                halfDual, halfWidth, halfHeight, 7);
        double* halfXGradient = new double[halfWidth * halfHeight];
        double* halfYGradient = new double[halfWidth * halfHeight];
        coefficientToXYGradient2D(halfImage, halfXGradient, halfYGradient,
                halfWidth, halfHeight);
        basicToCardinal2D(halfImage, halfImage, halfWidth, halfHeight, 3);
        pyramid.push(halfYGradient);
        pyramid.push(halfXGradient);
        pyramid.push(halfImage);
        pyramid.push(new Integer(halfHeight));
        pyramid.push(new Integer(halfWidth));
        fullDual = halfDual;
    }
} /* end buildImageAndGradientPyramid */

/*------------------------------------------------------------------*/
void TurboRegImage::buildImagePyramid (
) {
    int fullWidth;
    int fullHeight;
    double* fullDual = new double[width * height];
    int halfWidth = width;
    int halfHeight = height;
    if (1 < pyramidDepth) {
        cardinalToDual2D(image, fullDual, width, height, 3);
    }
    for (int depth = 1; ((depth < pyramidDepth));
            depth++) {
        fullWidth = halfWidth;
        fullHeight = halfHeight;
        halfWidth /= 2;
        halfHeight /= 2;
        double* halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
        double* halfImage = new double[halfWidth * halfHeight];
        dualToCardinal2D(halfDual, halfImage, halfWidth, halfHeight, 3);
        pyramid.push(halfImage);
        pyramid.push(new Integer(halfHeight));
        pyramid.push(new Integer(halfWidth));
        fullDual = halfDual;
    }
} /* end buildImagePyramid */

/*------------------------------------------------------------------*/
void TurboRegImage::cardinalToDual2D (
        double* cardinal,
        double* dual,
        int width,
        int height,
        int degree
) {
    basicToCardinal2D(getBasicFromCardinal2D(cardinal, width, height, degree),
            dual, width, height, 2 * degree + 1);
} /* end cardinalToDual2D */

/*------------------------------------------------------------------*/
void TurboRegImage::coefficientToGradient1D (
        double* c,
        int clen
) {
    double h[] = {0.0, 1.0 / 2.0};
    double* s = new double[clen]; //OK
    antiSymmetricFirMirrorOffBounds1D(h, sizeof(h), c, clen, s, clen);
    //System.arraycopy(s, 0, c, 0, slen);
    memcpy(c, s, slen * sizeof(double));
    delete s;
} /* end coefficientToGradient1D */

/*------------------------------------------------------------------*/
void TurboRegImage::coefficientToSamples1D (
        double* c,
        int clen
) {
    double h[] = {2.0 / 3.0, 1.0 / 6.0};
    double* s = new double[clen]; //OK
    symmetricFirMirrorOffBounds1D(h, COUNTOF(h), c, clen, s, clen);
    //System.arraycopy(s, 0, c, 0, slen);
    memcpy(c, s, slen * sizeof(double));
    delete s;
} /* end coefficientToSamples1D */

/*------------------------------------------------------------------*/
void TurboRegImage::coefficientToXYGradient2D (
        double* basic,
        double* xGradient,
        double* yGradient,
        int width,
        int height
) {
    double* hLine = new double[width];
    double* hData = new double[width];
    double* vLine = new double[height];
    
    for (int y = 0; ((y < height)); y++) {
        extractRow(basic, y, hLine, width);
        //System.arraycopy(hLine, 0, hData, 0, width);
        memcpy(hData, hLine, sizeof(double) * width);
        coefficientToGradient1D(hLine, width);
        
        coefficientToSamples1D(hData, width);
        putRow(xGradient, y, hLine, width);
        putRow(yGradient, y, hData, width);
    }
    
    for (int x = 0; ((x < width)); x++) {
        extractColumn(xGradient, width, x, vLine, height);
        coefficientToSamples1D(vLine, height);
        putColumn(xGradient, width, x, vLine, height);
        
        extractColumn(yGradient, width, x, vLine, height);
        coefficientToGradient1D(vLine, height);
        putColumn(yGradient, width, x, vLine, height);
    }
} /* end coefficientToXYGradient2D */

/*------------------------------------------------------------------*/
void TurboRegImage::dualToCardinal2D (
        double* dual,
        double* cardinal,
        int width,
        int height,
        int degree
) {
    basicToCardinal2D(
        getBasicFromCardinal2D(dual, width, height, 2 * degree + 1), 
            cardinal, width, height, degree);
} /* end dualToCardinal2D */

/*------------------------------------------------------------------*/
void TurboRegImage::extractColumn (
        double* array,
        int width,
        int x,
        double* column,
        int height
) {
    for (int i = 0; (i < height); i++) {
        column[i] = (double)array[x];
        x += width;
    }
} /* end extractColumn */

/*------------------------------------------------------------------*/
void TurboRegImage::extractRow (
        double* array,
        int y,
        double* row,
        int width
) {
    y *= width;
    for (int i = 0; (i < width); i++) {
        row[i] = (double)array[y++];
    }
} /* end extractRow */

/*------------------------------------------------------------------*/
double* TurboRegImage::getBasicFromCardinal2D (
) {
    double* basic = new double[width * height];
    double* hLine = new double[width];
    double* vLine = new double[height];
    
    for (int y = 0; (y < height); y++) {
        extractRow(image, y, hLine, width);
        samplesToInterpolationCoefficient1D(hLine, height, 3, 0.0);
        putRow(basic, y, hLine, width);
    }
    for (int x = 0; (x < width); x++) {
        extractColumn(basic, width, x, vLine, height);
        samplesToInterpolationCoefficient1D(vLine, height, 3, 0.0);
        putColumn(basic, width, x, vLine, height);
    }
    turboRegProgressBar.workloadDone(width + height);
    return(basic);
} /* end getBasicFromCardinal2D */

/*------------------------------------------------------------------*/
double* TurboRegImage::getBasicFromCardinal2D (
        double* cardinal,
        int width,
        int height,
        int degree
) {
    double* basic = new double[width * height];
    double* hLine = new double[width];
    double* vLine = new double[height];
    
    for (int y = 0; ((y < height)); y++) {
        extractRow(cardinal, y, hLine, width);
        samplesToInterpolationCoefficient1D(hLine, width, degree, 0.0);
        putRow(basic, y, hLine, width);
    }

    for (int x = 0; ((x < width)); x++) {
        extractColumn(basic, width, x, vLine, height);
        samplesToInterpolationCoefficient1D(vLine, height, degree, 0.0);
        putColumn(basic, width, x, vLine, height);
    }

    return(basic);
} /* end getBasicFromCardinal2D */

/*------------------------------------------------------------------*/
double* TurboRegImage::getHalfDual2D (
        double* fullDual,
        int fullWidth,
        int fullHeight
) {
    int halfWidth = fullWidth / 2;
    int halfHeight = fullHeight / 2;
    double* hLine = new double[fullWidth];
    double* hData = new double[halfWidth];
    double* vLine = new double[fullHeight];
    double* vData = new double[halfHeight];
    double* demiDual = new double[halfWidth * fullHeight];
    double* halfDual = new double[halfWidth * halfHeight];
    
    for (int y = 0; ((y < fullHeight)); y++) {
        extractRow(fullDual, y, hLine, width);
        reduceDual1D(hLine, width, hData, width);
        putRow(demiDual, y, hData, width);
    }

    for (int x = 0; ((x < halfWidth)); x++) {
        extractColumn(demiDual, halfWidth, x, vLine, height);
        reduceDual1D(vLine, height, vData, height);
        putColumn(halfDual, halfWidth, x, vData, height);
    }

    
    return(halfDual);
} /* end getHalfDual2D */

/*------------------------------------------------------------------*/
double TurboRegImage::getInitialAntiCausalCoefficientMirrorOffBounds (
        double* c,
        int clen,
        double z,
        double tolerance
) {
    return(z * c[clen - 1] / (z - 1.0));
} /* end getInitialAntiCausalCoefficientMirrorOffBounds */

/*------------------------------------------------------------------*/
double TurboRegImage::getInitialCausalCoefficientMirrorOffBounds (
        double* c,
        int clen,
        double z,
        double tolerance
) {
    double z1 = z;
    double zn = pow(z, clen);
    double sum = (1.0 + z) * (c[0] + zn * c[clen - 1]);
    int horizon = clen;
    if (0.0 < tolerance) {
        horizon = 2 + (int)(log(tolerance) / log((double)abs(z)));
        horizon = (horizon < clen) ? (horizon) : (clen);
    }
    zn = zn * zn;
    for (int n = 1; (n < (horizon - 1)); n++) {
        z1 = z1 * z;
        zn = zn / z;
        sum = sum + (z1 + zn) * c[n];
    }
    return(sum / (1.0 - pow(z, 2 * clen)));
} /* end getInitialCausalCoefficientMirrorOffBounds */

/*------------------------------------------------------------------*/
void TurboRegImage::imageToXYGradient2D (
) {
    double* hLine = new double[width];
    double* vLine = new double[height];
    xGradient = new double[width * height];
    yGradient = new double[width * height];
    int workload = width + height;
    turboRegProgressBar.addWorkload(workload);
    for (int y = 0; ((y < height)); y++) {
        extractRow(image, y, hLine, width);
        samplesToInterpolationCoefficient1D(hLine, width, 3, 0.0);
        coefficientToGradient1D(hLine, width);
        putRow(xGradient, y, hLine, width);
        turboRegProgressBar.stepProgressBar();
        workload--;
    }
    for (int x = 0; ((x < width)); x++) {
        extractColumn(image, width, x, vLine, height);
        samplesToInterpolationCoefficient1D(vLine, height, 3, 0.0);
        coefficientToGradient1D(vLine, height);
        putColumn(yGradient, width, x, vLine, height);
        turboRegProgressBar.stepProgressBar();
        workload--;
    }
    turboRegProgressBar.skipProgressBar(workload);
    turboRegProgressBar.workloadDone(width + height);
} /* end imageToXYGradient2D */

/*------------------------------------------------------------------*/
void TurboRegImage::putColumn (
        double* array,
        int width,
        int x,
        double* column,
        int columnlen
) {
    for (int i = 0; (i < columnlen); i++) {
        array[x] = (float)column[i];
        x += width;
    }
} /* end putColumn */

/*------------------------------------------------------------------*/
void TurboRegImage::putRow (
        double* array,
        int y,
        double* row,
        int rowlen
) {
    y *= rowlen;
    for (int i = 0; (i < rowlen); i++) {
        array[y++] = (float)row[i];
    }
} /* end putRow */

/*------------------------------------------------------------------*/
void TurboRegImage::reduceDual1D (
        double* c,
        int clen,
        double* s,
        int slen
) {
    double h[] = {6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0};
    if (2 <= slen) {
        s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]);
        for (int i = 2, j = 1; (j < (slen - 1)); i += 2, j++) {
            s[j] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1])
                    + h[2] * (c[i - 2] + c[i + 2]);
        }
        if (clen == (2 * slen)) {
            s[slen - 1] = h[0] * c[clen - 2]
                    + h[1] * (c[clen - 3] + c[clen - 1])
                    + h[2] * (c[clen - 4] + c[clen - 1]);
        }
        else {
            s[slen - 1] = h[0] * c[clen - 3]
                    + h[1] * (c[clen - 4] + c[clen - 2])
                    + h[2] * (c[clen - 5] + c[clen - 1]);
        }
    }
    else {
        switch (clen) {
            case 3: {
                s[0] = h[0] * c[0]
                        + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]);
                break;
            }
            case 2: {
                s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + 2.0 * h[2] * c[1];
                break;
            }
        }
    }
} /* end reduceDual1D */

/*------------------------------------------------------------------*/
void TurboRegImage::samplesToInterpolationCoefficient1D (
        double* c,
        int clen,
        int degree,
        double tolerance
) {
    double* z = new double[0];
    double lambda = 1.0;
    switch (degree) {
        case 3: {
            z = new double[1];
            z[0] = sqrt(3.0) - 2.0;
            break;
        }
        case 7: {
            z = new double[3];
            z[0] =
                    -0.5352804307964381655424037816816460718339231523426924148812;
            z[1] =
                    -0.122554615192326690515272264359357343605486549427295558490763;
            z[2] =
                    -0.0091486948096082769285930216516478534156925639545994482648003;
            break;
        }
    }
    if (clen == 1) {
        return;
    }

    for (int k = 0; (k < COUNTOF(z)); k++) {
        lambda *= (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
    }

    for (int n = 0; (n < clen); n++) {
        c[n] = c[n] * lambda;
    }

    for (int k = 0; (k < COUNTOF(z)); k++) {
        c[0] = getInitialCausalCoefficientMirrorOffBounds(c, clen, z[k], tolerance);
        for (int n = 1; (n < clen); n++) {
            c[n] = c[n] + z[k] * c[n - 1];
        }
        c[clen - 1] = getInitialAntiCausalCoefficientMirrorOffBounds(
                c, clen, z[k], tolerance);
        for (int n = clen - 2; (0 <= n); n--) {
            c[n] = z[k] * (c[n+1] - c[n]);
        }
    }
} /* end samplesToInterpolationCoefficient1D */

/*------------------------------------------------------------------*/
void TurboRegImage::symmetricFirMirrorOffBounds1D (
        double* h,
        int hlen,
        double* c,
        int clen,
        double* s,
        int slen
) {
    switch (hlen) {
        case 2: {
            if (2 <= clen) {
                s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]);
                for (int i = 1; (i < (slen - 1)); i++) {
                    s[i] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1]);
                }
                s[slen - 1] = h[0] * c[clen - 1]
                        + h[1] * (c[clen - 2] + c[clen - 1]);
            }
            else {
                s[0] = (h[0] + 2.0 * h[1]) * c[0];
            }
            break;
        }
        case 4: {
            if (6 <= clen) {
                s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2])
                        + h[3] * (c[2] + c[3]);
                s[1] = h[0] * c[1] + h[1] * (c[0] + c[2]) + h[2] * (c[0] + c[3])
                        + h[3] * (c[1] + c[4]);
                s[2] = h[0] * c[2] + h[1] * (c[1] + c[3]) + h[2] * (c[0] + c[4])
                        + h[3] * (c[0] + c[5]);
                for (int i = 3; (i < (slen - 3)); i++) {
                    s[i] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1])
                            + h[2] * (c[i - 2] + c[i + 2])
                            + h[3] * (c[i - 3] + c[i + 3]);
                }
                s[slen - 3] = h[0] * c[clen - 3]
                        + h[1] * (c[clen - 4] + c[clen - 2])
                        + h[2] * (c[clen - 5] + c[clen - 1])
                        + h[3] * (c[clen - 6] + c[clen - 1]);
                s[slen - 2] = h[0] * c[clen - 2]
                        + h[1] * (c[clen - 3] + c[clen - 1])
                        + h[2] * (c[clen - 4] + c[clen - 1])
                        + h[3] * (c[clen - 5] + c[clen - 2]);
                s[slen - 1] = h[0] * c[clen - 1]
                        + h[1] * (c[clen - 2] + c[clen - 1])
                        + h[2] * (c[clen - 3] + c[clen - 2])
                        + h[3] * (c[clen - 4] + c[clen - 3]);
            }
            else {
                switch (clen) {
                    case 5: {
                        s[0] = h[0] * c[0] + h[1] * (c[0] + c[1])
                                + h[2] * (c[1] + c[2]) + h[3] * (c[2] + c[3]);
                        s[1] = h[0] * c[1] + h[1] * (c[0] + c[2])
                                + h[2] * (c[0] + c[3]) + h[3] * (c[1] + c[4]);
                        s[2] = h[0] * c[2] + h[1] * (c[1] + c[3])
                                + (h[2] + h[3]) * (c[0] + c[4]);
                        s[3] = h[0] * c[3] + h[1] * (c[2] + c[4])
                                + h[2] * (c[1] + c[4]) + h[3] * (c[0] + c[3]);
                        s[4] = h[0] * c[4] + h[1] * (c[3] + c[4])
                                + h[2] * (c[2] + c[3]) + h[3] * (c[1] + c[2]);
                        break;
                    }
                    case 4: {
                        s[0] = h[0] * c[0] + h[1] * (c[0] + c[1])
                                + h[2] * (c[1] + c[2]) + h[3] * (c[2] + c[3]);
                        s[1] = h[0] * c[1] + h[1] * (c[0] + c[2])
                                + h[2] * (c[0] + c[3]) + h[3] * (c[1] + c[3]);
                        s[2] = h[0] * c[2] + h[1] * (c[1] + c[3])
                                + h[2] * (c[0] + c[3]) + h[3] * (c[0] + c[2]);
                        s[3] = h[0] * c[3] + h[1] * (c[2] + c[3])
                                + h[2] * (c[1] + c[2]) + h[3] * (c[0] + c[1]);
                        break;
                    }
                    case 3: {
                        s[0] = h[0] * c[0] + h[1] * (c[0] + c[1])
                                + h[2] * (c[1] + c[2]) + 2.0 * h[3] * c[2];
                        s[1] = h[0] * c[1] + (h[1] + h[2]) * (c[0] + c[2])
                                + 2.0 * h[3] * c[1];
                        s[2] = h[0] * c[2] + h[1] * (c[1] + c[2])
                                + h[2] * (c[0] + c[1]) + 2.0 * h[3] * c[0];
                        break;
                    }
                    case 2: {
                        s[0] = (h[0] + h[1] + h[3]) * c[0]
                                + (h[1] + 2.0 * h[2] + h[3]) * c[1];
                        s[1] = (h[0] + h[1] + h[3]) * c[1]
                                + (h[1] + 2.0 * h[2] + h[3]) * c[0];
                        break;
                    }
                    case 1: {
                        s[0] = (h[0] + 2.0 * (h[1] + h[2] + h[3])) * c[0];
                        break;
                    }
                }
            }
            break;
        }
    }
} /* end symmetricFirMirrorOffBounds1D */

