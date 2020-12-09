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

#include <cstring>
#include <vector>
#include <cmath>

#ifdef _MSC_VER
template <typename T, size_t N>
T* begin(T(&arr)[N]) { return &arr[0]; }
template <typename T, size_t N>
T* end(T(&arr)[N]) { return &arr[0] + N; }
#endif


//#define COUNTOF(x) (sizeof(x) / sizeof(*x))

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
	this->transformation = transformation;
    this->isTarget = isTarget;
    this->width = width;
    this->height = height;
    int k = 0;

    image.resize(width * height);
    for (int y = 0; (y < height); y++) {
        for (int x = 0; (x < width); x++, k++) {
            image[k] = (double)img[k];
        }
    }



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

TurboRegImage::TurboRegImage (
    double *img, int width, int height, bool isTarget
) {
    this->isTarget = isTarget;
    this->width = width;
    this->height = height;
    int k = 0;

    image.resize(width * height);
    for (int y = 0; (y < height); y++) {
        for (int x = 0; (x < width); x++, k++) {
            image[k] = (double)img[k];
        }
    }
}

void TurboRegImage::antiSymmetricFirMirrorOffBounds1D(
        std::vector<double> &h,
        std::vector<double> &c,
        std::vector<double> &s
) {
    if (2 <= h.size()) {
        s[0] = h[1] * (c[1] - c[0]);
        for (unsigned int i = 1; (i < (s.size() - 1)); i++) {
            s[i] = h[1] * (c[i + 1] - c[i - 1]);
        }
        s[s.size() - 1] = h[1] * (c[c.size() - 1] - c[c.size() - 2]);
    }
    else {
        s[0] = 0.0;
    }
} /* end antiSymmetricFirMirrorOffBounds1D */

/*------------------------------------------------------------------*/
void TurboRegImage::basicToCardinal2D (
        const std::vector<double> &basic,
        std::vector<double> &cardinal,
        int width,
        int height,
        int degree
) {
    std::vector<double> hLine(width);
    std::vector<double> vLine(height);
    std::vector<double> hData(width);
    std::vector<double> vData(height);

    std::vector<double> h;

    switch (degree) {
        case 3: {
            //h = new double[2];
            h.resize(2);
            h[0] = 2.0 / 3.0;
            h[1] = 1.0 / 6.0;
            break;
        }
        case 7: {
            //h = new double[4];
            h.resize(4);
            h[0] = 151.0 / 315.0;
            h[1] = 397.0 / 1680.0;
            h[2] = 1.0 / 42.0;
            h[3] = 1.0 / 5040.0;
            break;
        }
        default: {
            h.resize(1);
            //h = new double[1];
            h[0] = 1.0;
        }
    }

    for (int y = 0; (y < height); y++) {
        extractRow(basic, y, hLine);
        symmetricFirMirrorOffBounds1D(h, hLine, hData);
        putRow(cardinal, y, hData);
    }

    for (int x = 0; (x < width); x++) {
        extractColumn(cardinal, width, x, vLine);
        symmetricFirMirrorOffBounds1D(h, vLine, vData);
        putColumn(cardinal, width, x, vData);
    }

} /* end basicToCardinal2D */

/*------------------------------------------------------------------*/
void TurboRegImage::buildCoefficientPyramid (
) {
    int fullWidth;
    int fullHeight;
    std::vector<double> fullDual(width * height);
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
        std::vector<double> halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);

        ImageStackItem stackItem(halfWidth, halfHeight, false);

        //std::vector<double> halfCoefficient =
        getBasicFromCardinal2D(
                halfDual, halfWidth, halfHeight, 7, stackItem.halfImg);
        pyramid.push(stackItem);
        /*pyramid.push(halfCoefficient);
        pyramid.push(new Integer(halfHeight));
        pyramid.push(new Integer(halfWidth));*/

        fullDual = halfDual;
    }
} /* end buildCoefficientPyramid */

/*------------------------------------------------------------------*/
void TurboRegImage::buildImageAndGradientPyramid (
) {
    int fullWidth;
    int fullHeight;
    std::vector<double> fullDual(width * height);
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

        ImageStackItem stackItem(halfWidth, halfHeight, true);

        std::vector<double> halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
        //std::vector<double> halfImage =
        getBasicFromCardinal2D(
                halfDual, halfWidth, halfHeight, 7, stackItem.halfImg);
        //std::vector<double> halfXGradient(halfWidth * halfHeight);
        //std::vector<double> halfYGradient(halfWidth * halfHeight);
        coefficientToXYGradient2D(stackItem.halfImg, stackItem.xGradient, stackItem.yGradient,
                halfWidth, halfHeight);
        basicToCardinal2D(stackItem.halfImg, stackItem.halfImg, halfWidth, halfHeight, 3);

        pyramid.push(stackItem);

        /*pyramid.push(halfYGradient);
        pyramid.push(halfXGradient);
        pyramid.push(halfImage);
        pyramid.push(new Integer(halfHeight));
        pyramid.push(new Integer(halfWidth));*/
        fullDual = halfDual;
    }
} /* end buildImageAndGradientPyramid */

/*------------------------------------------------------------------*/
void TurboRegImage::buildImagePyramid (
) {
    int fullWidth;
    int fullHeight;
    std::vector<double> fullDual(width * height);
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

        ImageStackItem stackItem(halfWidth, halfHeight, true);


        std::vector<double> halfDual = getHalfDual2D(fullDual, fullWidth, fullHeight);
        //std::vector<double> halfImage(halfWidth * halfHeight);
        dualToCardinal2D(halfDual, stackItem.halfImg, halfWidth, halfHeight, 3);
        /*pyramid.push(halfImage);
        pyramid.push(new Integer(halfHeight));
        pyramid.push(new Integer(halfWidth));*/

        pyramid.push(stackItem);

        fullDual = halfDual;
    }
} /* end buildImagePyramid */

/*------------------------------------------------------------------*/
void TurboRegImage::cardinalToDual2D (
        std::vector<double> &cardinal,
        std::vector<double> &dual,
        int width,
        int height,
        int degree
) {
    std::vector<double> basic(width * height);
    basicToCardinal2D(getBasicFromCardinal2D(cardinal, width, height, degree, basic),
            dual, width, height, 2 * degree + 1);
} /* end cardinalToDual2D */

/*------------------------------------------------------------------*/
void TurboRegImage::coefficientToGradient1D (
        std::vector<double> &c
) {
#ifdef _MSC_VER
	double hh[] = { 0.0, 1.0 / 2.0 };
	std::vector<double> h(begin(hh), end(hh));
#else
	std::vector<double> h = { 0.0, 1.0 / 2.0 };
#endif

    std::vector<double> s(c.size()); //OK
    antiSymmetricFirMirrorOffBounds1D(h, c, s);
    //System.arraycopy(s, 0, c, 0, s.size());
    c = s;
} /* end coefficientToGradient1D */

/*------------------------------------------------------------------*/
void TurboRegImage::coefficientToSamples1D (
        std::vector<double> &c
) {
#ifdef _MSC_VER
	double hh[] = { 2.0 / 3.0, 1.0 / 6.0 };
	std::vector<double> h(begin(hh), end(hh));
#else
	std::vector<double> h = { 2.0 / 3.0, 1.0 / 6.0 };
#endif

    std::vector<double> s(c.size()); //OK
    symmetricFirMirrorOffBounds1D(h, c, s);
    //System.arraycopy(s, 0, c, 0, s.size());
    c = s;
} /* end coefficientToSamples1D */

/*------------------------------------------------------------------*/
void TurboRegImage::coefficientToXYGradient2D (
        std::vector<double> &basic,
        std::vector<double> &xGradient,
        std::vector<double> &yGradient,
        int width,
        int height
) {
    std::vector<double> hLine(width);
    std::vector<double> hData(width);
    std::vector<double> vLine(height);

    for (int y = 0; ((y < height)); y++) {
        extractRow(basic, y, hLine);
        //System.arraycopy(hLine, 0, hData, 0, width);
        //memcpy(hData, hLine, sizeof(double) * width);
        hData = hLine;
        coefficientToGradient1D(hLine);

        coefficientToSamples1D(hData);
        putRow(xGradient, y, hLine);
        putRow(yGradient, y, hData);
    }

    for (int x = 0; ((x < width)); x++) {
        extractColumn(xGradient, width, x, vLine);
        coefficientToSamples1D(vLine);
        putColumn(xGradient, width, x, vLine);

        extractColumn(yGradient, width, x, vLine);
        coefficientToGradient1D(vLine);
        putColumn(yGradient, width, x, vLine);
    }
} /* end coefficientToXYGradient2D */

/*------------------------------------------------------------------*/
void TurboRegImage::dualToCardinal2D (
        std::vector<double> &dual,
        std::vector<double> &cardinal,
        int width,
        int height,
        int degree
) {
    std::vector<double> basic(width * height);
    basicToCardinal2D(
        getBasicFromCardinal2D(dual, width, height, 2 * degree + 1, basic),
            cardinal, width, height, degree);
} /* end dualToCardinal2D */

/*------------------------------------------------------------------*/
void TurboRegImage::extractColumn (
        std::vector<double> &array,
        int width,
        int x,
        std::vector<double> &column
) {
    for (unsigned int i = 0; (i < column.size()); i++) {
        column[i] = (double)array[x];
        x += width;
    }
} /* end extractColumn */

/*------------------------------------------------------------------*/
void TurboRegImage::extractRow (
        const std::vector<double> &array,
        int y,
        std::vector<double> &row
) {
    y *= (int)row.size();
    for (unsigned int i = 0; (i < row.size()); i++) {
        row[i] = (double)array[y++];
    }
} /* end extractRow */

/*------------------------------------------------------------------*/
std::vector<double> TurboRegImage::getBasicFromCardinal2D (
) {
    std::vector<double> basic(width * height);
    std::vector<double> hLine(width);
    std::vector<double> vLine(height);

    for (int y = 0; (y < height); y++) {
        extractRow(image, y, hLine);
        samplesToInterpolationCoefficient1D(hLine, 3, 0.0);
        putRow(basic, y, hLine);
    }
    for (int x = 0; (x < width); x++) {
        extractColumn(basic, width, x, vLine);
        samplesToInterpolationCoefficient1D(vLine, 3, 0.0);
        putColumn(basic, width, x, vLine);
    }

    return(basic);
} /* end getBasicFromCardinal2D */

/*------------------------------------------------------------------*/
std::vector<double> TurboRegImage::getBasicFromCardinal2D (
        std::vector<double> &cardinal,
        int width,
        int height,
        int degree,
        std::vector<double> &basic
) {
    //std::vector<double> basic(width * height);
    std::vector<double> hLine(width);
    std::vector<double> vLine(height);

    for (int y = 0; ((y < height)); y++) {
        extractRow(cardinal, y, hLine);
        samplesToInterpolationCoefficient1D(hLine, degree, 0.0);
        putRow(basic, y, hLine);
    }

    for (int x = 0; ((x < width)); x++) {
        extractColumn(basic, width, x, vLine);
        samplesToInterpolationCoefficient1D(vLine, degree, 0.0);
        putColumn(basic, width, x, vLine);
    }

    return(basic);
} /* end getBasicFromCardinal2D */

/*------------------------------------------------------------------*/
std::vector<double> TurboRegImage::getHalfDual2D (
        std::vector<double> &fullDual,
        int fullWidth,
        int fullHeight
) {
    int halfWidth = fullWidth / 2;
    int halfHeight = fullHeight / 2;
    std::vector<double> hLine(fullWidth);
    std::vector<double> hData(halfWidth);
    std::vector<double> vLine(fullHeight);
    std::vector<double> vData(halfHeight);
    std::vector<double> demiDual(halfWidth * fullHeight);
    std::vector<double> halfDual(halfWidth * halfHeight);

    for (int y = 0; ((y < fullHeight)); y++) {
        extractRow(fullDual, y, hLine);
        reduceDual1D(hLine, hData);
        putRow(demiDual, y, hData);
    }

    for (int x = 0; ((x < halfWidth)); x++) {
        extractColumn(demiDual, halfWidth, x, vLine);
        reduceDual1D(vLine, vData);
        putColumn(halfDual, halfWidth, x, vData);
    }


    return(halfDual);
} /* end getHalfDual2D */

/*------------------------------------------------------------------*/
double TurboRegImage::getInitialAntiCausalCoefficientMirrorOffBounds (
        std::vector<double> &c,
        double z,
        double tolerance
) {
    return(z * c[c.size() - 1] / (z - 1.0));
} /* end getInitialAntiCausalCoefficientMirrorOffBounds */

/*------------------------------------------------------------------*/
double TurboRegImage::getInitialCausalCoefficientMirrorOffBounds (
        std::vector<double> &c,
        double z,
        double tolerance
) {
    double z1 = z;
    double zn = pow(z, (double)c.size());
    double sum = (1.0 + z) * (c[0] + zn * c[c.size() - 1]);
    int horizon = (int)c.size();
    if (0.0 < tolerance) {
        horizon = 2 + (int)(log(tolerance) / log((double)std::abs(z)));
        horizon = (horizon < (int)c.size()) ? (horizon) : ((int)c.size());
    }
    zn = zn * zn;
    for (int n = 1; (n < (horizon - 1)); n++) {
        z1 = z1 * z;
        zn = zn / z;
        sum = sum + (z1 + zn) * c[n];
    }
    return(sum / (1.0 - pow(z, (double)(2 * c.size()))));
} /* end getInitialCausalCoefficientMirrorOffBounds */

/*------------------------------------------------------------------*/
void TurboRegImage::imageToXYGradient2D (
) {
    std::vector<double> hLine(width);
    std::vector<double> vLine(height);
    xGradient.resize(width * height);
    yGradient.resize(width * height);

    for (int y = 0; ((y < height)); y++) {
        extractRow(image, y, hLine);
        samplesToInterpolationCoefficient1D(hLine, 3, 0.0);
        coefficientToGradient1D(hLine);
        putRow(xGradient, y, hLine);
    }
    for (int x = 0; ((x < width)); x++) {
        extractColumn(image, width, x, vLine);
        samplesToInterpolationCoefficient1D(vLine, 3, 0.0);
        coefficientToGradient1D(vLine);
        putColumn(yGradient, width, x, vLine);
    }
} /* end imageToXYGradient2D */

/*------------------------------------------------------------------*/
void TurboRegImage::putColumn (
        std::vector<double> &array,
        int width,
        int x,
        std::vector<double> &column
) {
    for (unsigned int i = 0; (i < column.size()); i++) {
        array[x] = (float)column[i];
        x += width;
    }
} /* end putColumn */

/*------------------------------------------------------------------*/
void TurboRegImage::putRow (
        std::vector<double> &array,
        int y,
        std::vector<double> &row
) {
    y *= row.size();
    for (unsigned int i = 0; (i < row.size()); i++) {
        array[y++] = (float)row[i];
    }
} /* end putRow */

/*------------------------------------------------------------------*/
void TurboRegImage::reduceDual1D (
        std::vector<double> &c,
        std::vector<double> &s
) {
    double h[] = {6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0};
    if (2 <= s.size()) {
        s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2]);
        for (unsigned int i = 2, j = 1; (j < (s.size() - 1)); i += 2, j++) {
            s[j] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1])
                    + h[2] * (c[i - 2] + c[i + 2]);
        }
        if (c.size() == (2 * s.size())) {
            s[s.size() - 1] = h[0] * c[c.size() - 2]
                    + h[1] * (c[c.size() - 3] + c[c.size() - 1])
                    + h[2] * (c[c.size() - 4] + c[c.size() - 1]);
        }
        else {
            s[s.size() - 1] = h[0] * c[c.size() - 3]
                    + h[1] * (c[c.size() - 4] + c[c.size() - 2])
                    + h[2] * (c[c.size() - 5] + c[c.size() - 1]);
        }
    }
    else {
        switch (c.size()) {
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
        std::vector<double> &c,
        int degree,
        double tolerance
) {
    std::vector<double> z;
    double lambda = 1.0;
    switch (degree) {
        case 3: {
            //z = new double[1];
            z.resize(1);

            z[0] = sqrt(3.0) - 2.0;
            break;
        }
        case 7: {
            //z = new double[3];
            z.resize(3);
            z[0] =
                    -0.5352804307964381655424037816816460718339231523426924148812;
            z[1] =
                    -0.122554615192326690515272264359357343605486549427295558490763;
            z[2] =
                    -0.0091486948096082769285930216516478534156925639545994482648003;
            break;
        }
    }
    if (c.size() == 1) {
        return;
    }

    for (unsigned int k = 0; (k < z.size()); k++) {
        lambda *= (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
    }

    for (unsigned int n = 0; (n < c.size()); n++) {
        c[n] = c[n] * lambda;
    }

    for (unsigned int k = 0; (k < z.size()); k++) {
        c[0] = getInitialCausalCoefficientMirrorOffBounds(c, z[k], tolerance);
        for (unsigned int n = 1; (n < c.size()); n++) {
            c[n] = c[n] + z[k] * c[n - 1];
        }
        c[c.size() - 1] = getInitialAntiCausalCoefficientMirrorOffBounds(
                c, z[k], tolerance);
        for (int n = c.size() - 2; (0 <= n); n--) {
            c[n] = z[k] * (c[n+1] - c[n]);
        }
    }
} /* end samplesToInterpolationCoefficient1D */

/*------------------------------------------------------------------*/
void TurboRegImage::symmetricFirMirrorOffBounds1D (
        std::vector<double> &h,
        std::vector<double> &c,
        std::vector<double> &s
) {
    switch (h.size()) {
        case 2: {
            if (2 <= c.size()) {
                s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]);
                for (unsigned int i = 1; (i < (s.size() - 1)); i++) {
                    s[i] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1]);
                }
                s[s.size() - 1] = h[0] * c[c.size() - 1]
                        + h[1] * (c[c.size() - 2] + c[c.size() - 1]);
            }
            else {
                s[0] = (h[0] + 2.0 * h[1]) * c[0];
            }
            break;
        }
        case 4: {
            if (6 <= c.size()) {
                s[0] = h[0] * c[0] + h[1] * (c[0] + c[1]) + h[2] * (c[1] + c[2])
                        + h[3] * (c[2] + c[3]);
                s[1] = h[0] * c[1] + h[1] * (c[0] + c[2]) + h[2] * (c[0] + c[3])
                        + h[3] * (c[1] + c[4]);
                s[2] = h[0] * c[2] + h[1] * (c[1] + c[3]) + h[2] * (c[0] + c[4])
                        + h[3] * (c[0] + c[5]);
                for (unsigned int i = 3; (i < (s.size() - 3)); i++) {
                    s[i] = h[0] * c[i] + h[1] * (c[i - 1] + c[i + 1])
                            + h[2] * (c[i - 2] + c[i + 2])
                            + h[3] * (c[i - 3] + c[i + 3]);
                }
                s[s.size() - 3] = h[0] * c[c.size() - 3]
                        + h[1] * (c[c.size() - 4] + c[c.size() - 2])
                        + h[2] * (c[c.size() - 5] + c[c.size() - 1])
                        + h[3] * (c[c.size() - 6] + c[c.size() - 1]);
                s[s.size() - 2] = h[0] * c[c.size() - 2]
                        + h[1] * (c[c.size() - 3] + c[c.size() - 1])
                        + h[2] * (c[c.size() - 4] + c[c.size() - 1])
                        + h[3] * (c[c.size() - 5] + c[c.size() - 2]);
                s[s.size() - 1] = h[0] * c[c.size() - 1]
                        + h[1] * (c[c.size() - 2] + c[c.size() - 1])
                        + h[2] * (c[c.size() - 3] + c[c.size() - 2])
                        + h[3] * (c[c.size() - 4] + c[c.size() - 3]);
            }
            else {
                switch (c.size()) {
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
