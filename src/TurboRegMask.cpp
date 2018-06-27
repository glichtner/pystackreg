#include "TurboRegMask.h"
#include <cmath>

/*====================================================================
|	turboRegMask
\===================================================================*/

/*********************************************************************
 This class is responsible for the mask preprocessing that takes
place concurrently with user-interface events. It contains methods
to compute the mask pyramids.
********************************************************************/



/*....................................................................
Runnable methods
....................................................................*/
/*********************************************************************
 Start the mask precomputations, which are interruptible.
    ********************************************************************/
void TurboRegMask::init (
) {
    this->buildPyramid();
} /* end run */


/*....................................................................
constructors
....................................................................*/
/*********************************************************************
 Converts the pixel array of the incoming <code>ImagePlus</code>
    object into a local <code>boolean</code> array.
    @param imp <code>ImagePlus</code> object to preprocess.
    ********************************************************************/
TurboRegMask::TurboRegMask (
        TurboRegImage &img
) {
    this->width = img.getWidth();
    this->height = img.getHeight();
    int k = 0;
    
    this->mask.resize(width * height);

    for (int y = 0; (y < height); y++) {
        for (int x = 0; (x < width); x++, k++) {
            mask[k] = img.getImage()[k];
        }
    }
}

TurboRegMask::TurboRegMask (
        matrix<double> &imp, int width, int height
) {
    this->width = width;   //imp.getWidth();
    this->height = height; //imp.getHeight();
    int k = 0;
    
    this->mask.resize(width * height);

    for (int y = 0; (y < height); y++) {
        for (int x = 0; (x < width); x++, k++) {
            mask[k] = (double)imp(x,y); //imp[k];
        }
    }

    /*if (imp.getType() == ImagePlus.GRAY8) {
        byte[] pixels = (byte[])imp.getProcessor().getPixels();
        for (int y = 0; (y < height); y++) {
            for (int x = 0; (x < width); x++, k++) {
                mask[k] = (float)pixels[k];
            }
        }
    }
    else if (imp.getType() == ImagePlus.GRAY16) {
        short[] pixels = (short[])imp.getProcessor().getPixels();
        for (int y = 0; (y < height); y++) {
            for (int x = 0; (x < width); x++, k++) {
                mask[k] = (float)pixels[k];
            }
        }
    }
    else if (imp.getType() == ImagePlus.GRAY32) {
        double* pixels = (double*)imp.getProcessor().getPixels();
        for (int y = 0; (y < height); y++) {
            for (int x = 0; (x < width); x++, k++) {
                mask[k] = pixels[k];
            }
        }
    }
    */


} /* end turboRegMask */

/*....................................................................
methods
....................................................................*/
/*********************************************************************
 Set to <code>true</code> every pixel of the full-size mask.
    ********************************************************************/
void TurboRegMask::clearMask (
) {
    int k = 0;
    for (int y = 0; (y < this->height); y++) {
        for (int x = 0; (x < this->width); x++) {
            this->mask[k++] = 1.0F;
        }
    }
} /* end clearMask */

/*********************************************************************
 Return the full-size mask array.
    ********************************************************************/
std::vector<double> &TurboRegMask::getMask (
) {
    return(this->mask);
} /* end getMask */

/*********************************************************************
 Return the pyramid as a <code>Stack</code> object. A single pyramid
    level consists of
    <p>
    <table border="1">
    <tr><th><code>isTarget</code></th>
    <th>ML*</th>
    <th>ML</th></tr>
    <tr><td>true</td>
    <td>mask samples</td>
    <td>mask samples</td></tr>
    <tr><td>false</td>
    <td>mask samples</td>
    <td>mask samples</td></tr>
    </table>
    @see turboRegImage#getPyramid()
    ********************************************************************/
std::stack<MaskStackItem> TurboRegMask::getPyramid (
) {
    return(this->pyramid);
} /* end getPyramid */


/*********************************************************************
 Set the depth up to which the pyramids should be computed.
    @see turboRegMask#getPyramid()
    ********************************************************************/
void TurboRegMask::setPyramidDepth (
        int pyramidDepth
) {
    this->pyramidDepth = pyramidDepth;
} /* end setPyramidDepth */

/*....................................................................
methods
....................................................................*/
/*------------------------------------------------------------------*/
void TurboRegMask::buildPyramid (
) {
    int fullWidth;
    int fullHeight;
    std::vector<double> &fullMask = mask;
    int halfWidth = width;
    int halfHeight = height;
    for (int depth = 1; ((depth < pyramidDepth)); depth++) {
        fullWidth = halfWidth;
        fullHeight = halfHeight;
        halfWidth /= 2;
        halfHeight /= 2;

        MaskStackItem stackItem(halfWidth * halfHeight);

        //std::vector<double> halfMask = getHalfMask2D(fullMask, fullWidth, fullHeight);
        getHalfMask2D(fullMask, fullWidth, fullHeight, stackItem.halfMask);

        pyramid.push(stackItem);
        fullMask = stackItem.halfMask;
    }
} /* end buildPyramid */

/*------------------------------------------------------------------*/
std::vector<double> TurboRegMask::getHalfMask2D (
        std::vector<double> &fullMask,
        int fullWidth,
        int fullHeight,
        std::vector<double> &halfMask
) {
    int halfWidth = fullWidth / 2;
    int halfHeight = fullHeight / 2;
    bool oddWidth = ((2 * halfWidth) != fullWidth);

    //std::vector<double> halfMask(halfWidth * halfHeight);
    int k = 0;
    for (int y = 0; y < halfHeight; y++) {
        for (int x = 0; (x < halfWidth); x++) {
            halfMask[k++] = 0.0F;
        }
    }

    k = 0;
    int n = 0;
    for (int y = 0; y < (halfHeight - 1); y++) {
        for (int x = 0; (x < (halfWidth - 1)); x++) {
            halfMask[k] += std::abs(fullMask[n++]);
            halfMask[k] += std::abs(fullMask[n]);
            halfMask[++k] += std::abs(fullMask[n++]);
        }
        halfMask[k] += std::abs(fullMask[n++]);
        halfMask[k++] += std::abs(fullMask[n++]);
        
        if (oddWidth) {
            n++;
        }
        for (int x = 0; (x < (halfWidth - 1)); x++) {
            halfMask[k - halfWidth] += std::abs(fullMask[n]);
            halfMask[k] += std::abs(fullMask[n++]);
            halfMask[k - halfWidth] += std::abs(fullMask[n]);
            halfMask[k - halfWidth + 1] += std::abs(fullMask[n]);
            halfMask[k] += std::abs(fullMask[n]);
            halfMask[++k] += std::abs(fullMask[n++]);
        }
        halfMask[k - halfWidth] += std::abs(fullMask[n]);
        halfMask[k] += std::abs(fullMask[n++]);
        halfMask[k - halfWidth] += std::abs(fullMask[n]);
        halfMask[k++] += std::abs(fullMask[n++]);
        
        if (oddWidth) {
            n++;
        }
        k -= halfWidth;
        
    }

    for (int x = 0; (x < (halfWidth - 1)); x++) {
        halfMask[k] += std::abs(fullMask[n++]);
        halfMask[k] += std::abs(fullMask[n]);
        halfMask[++k] += std::abs(fullMask[n++]);
    }

    halfMask[k] += std::abs(fullMask[n++]);
    halfMask[k++] += std::abs(fullMask[n++]);
    
    if (oddWidth) {
        n++;
    }
    
    k -= halfWidth;

    for (int x = 0; (x < (halfWidth - 1)); x++) {
        halfMask[k] += std::abs(fullMask[n++]);
        halfMask[k] += std::abs(fullMask[n]);
        halfMask[++k] += std::abs(fullMask[n++]);
    }
    halfMask[k] += std::abs(fullMask[n++]);
    halfMask[k] += std::abs(fullMask[n]);
    return(halfMask);
} /* end getHalfMask2D */

