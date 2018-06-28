#include <stack>
#include <vector>
#include "TurboReg.h"

#ifndef TURBOREGIMAGE_H_
#define TURBOREGIMAGE_H_

class ImageStackItem {
public:
    std::vector<double> halfImg;
    std::vector<double> xGradient;
    std::vector<double> yGradient;
    int halfWidth;
    int halfHeight;
    
    

    ImageStackItem(int halfWidth, int halfHeight, bool gradient)
    : halfWidth{halfWidth}, halfHeight{halfHeight}
     {
        halfImg.resize(halfHeight * halfWidth);

        if (gradient) {
            xGradient.resize(halfHeight * halfWidth);
            yGradient.resize(halfHeight * halfWidth);
        }
    }
};

/*====================================================================
|	turboRegImage
\===================================================================*/

/*********************************************************************
 This class is responsible for the image preprocessing that takes
 place concurrently with user-interface events. It contains methods
 to compute B-spline coefficients and their pyramids, image pyramids,
 gradients, and gradient pyramids.
 ********************************************************************/
class TurboRegImage { /* class turboRegImage */

public:
    TurboRegImage (double *img, int width, int height, Transformation transformation, bool isTarget);
    TurboRegImage (double *img, int width, int height, bool isTarget);
    
    std::vector<double> &getCoefficient () {
        return(this->coefficient);
    };

    int getHeight () {
        return(this->height);
    };

    std::vector<double> &getImage () {
        return(this->image);
    };

    std::stack<ImageStackItem> &getPyramid (
    ) {
        return (this->pyramid);
    };

    int getPyramidDepth (
    ) {
        return(this->pyramidDepth);
    };

    int getWidth (
    ) {
        return(this->width);
    } ;

    std::vector<double> &getXGradient (
    ) {
        return(this->xGradient);
    };

    std::vector<double> &getYGradient (
    ) {
        return(this->yGradient);
    };

    void setPyramidDepth ( int pyramidDepth) {
        this->pyramidDepth = pyramidDepth;
    }

    void setTransformation (Transformation transformation) {
        this->transformation = transformation;
    }

    void init();


private:
    std::stack<ImageStackItem> pyramid;
    
    std::vector<double> image;
    std::vector<double> coefficient;
    std::vector<double> xGradient;
    std::vector<double> yGradient;
    int width;
    int height;
    int pyramidDepth;
    Transformation transformation;
    bool isTarget;



    void antiSymmetricFirMirrorOffBounds1D (
        std::vector<double> &h,
        std::vector<double> &c,
        std::vector<double> &s
    );

    void basicToCardinal2D (
            const std::vector<double> &basic,
            std::vector<double> &cardinal,
            int width,
            int height,
            int degree
    );
    void buildCoefficientPyramid ();
    void buildImageAndGradientPyramid ();
    void buildImagePyramid ();
    void cardinalToDual2D (
            std::vector<double> &cardinal,
            std::vector<double> &dual,
            int width,
            int height,
            int degree
    );

    /*------------------------------------------------------------------*/
    void coefficientToGradient1D (
            std::vector<double> &c
    );

    /*------------------------------------------------------------------*/
    void coefficientToSamples1D (
            std::vector<double> &c
    ) ;

    /*------------------------------------------------------------------*/
    void coefficientToXYGradient2D (
            std::vector<double> &basic,
            std::vector<double> &xGradient,
            std::vector<double> &yGradient,
            int width,
            int height
    );

    /*------------------------------------------------------------------*/
    void dualToCardinal2D (
            std::vector<double> &dual,
            std::vector<double> &cardinal,
            int width,
            int height,
            int degree
    );

    /*------------------------------------------------------------------*/
    void extractColumn (
            std::vector<double> &array,
            int width,
            int x,
            std::vector<double> &column
    );

    /*------------------------------------------------------------------*/
    void extractRow (
        const std::vector<double> &array,
        int y,
        std::vector<double> &row
    );

    /*------------------------------------------------------------------*/
    std::vector<double> getBasicFromCardinal2D (
    );

    /*------------------------------------------------------------------*/
    std::vector<double> getBasicFromCardinal2D (
            std::vector<double> &cardinal,
            int width,
            int height,
            int degree,
            std::vector<double> &img
    );
    /*------------------------------------------------------------------*/
    std::vector<double> getHalfDual2D (
            std::vector<double> &fullDual,
            int fullWidth,
            int fullHeight
    );

    /*------------------------------------------------------------------*/
    double getInitialAntiCausalCoefficientMirrorOffBounds (
            std::vector<double> &c,
            double z,
            double tolerance
    );

    /*------------------------------------------------------------------*/
    double getInitialCausalCoefficientMirrorOffBounds (
            std::vector<double> &c,
            double z,
            double tolerance
    );

    void imageToXYGradient2D (
    );

    /*------------------------------------------------------------------*/
    void putColumn (
            std::vector<double> &array,
            int width,
            int x,
            std::vector<double> &column
    );
    /*------------------------------------------------------------------*/
    void putRow (
        std::vector<double> &array,
        int y,
        std::vector<double> &row
    );
    /*------------------------------------------------------------------*/
    void reduceDual1D (
        std::vector<double> &c,
        std::vector<double> &s
    );
    
    void samplesToInterpolationCoefficient1D (
            std::vector<double> &c,
            int degree,
            double tolerance
    );

    void symmetricFirMirrorOffBounds1D (
        std::vector<double> &h,
        std::vector<double> &c,
        std::vector<double> &s
    );

};

#endif
