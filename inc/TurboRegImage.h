#include <stack>
#include "TurboReg.h"

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
    
    double* getCoefficient () {
        return(this->coefficient);
    };

    int getHeight () {
        return(this->height);
    };

    double* getImage () {
        return(this->image);
    };

    Stack<Object> getPyramid (
    ) {
        return(this->pyramid);
    };

    int getPyramidDepth (
    ) {
        return(this->pyramidDepth);
    };

    int getWidth (
    ) {
        return(this->width);
    } ;

    double* getXGradient (
    ) {
        return(this->xGradient);
    };

    double* getYGradient (
    ) {
        return(this->yGradient);
    };

    void setPyramidDepth ( int pyramidDepth) {
        this->pyramidDepth = pyramidDepth;
    }

    void setTransformation (Transformation transformation) {
        this->transformation = transformation;
    }


private:
    std::stack<int> pyramid;
    
    double* image;
    double* coefficient;
    double* xGradient;
    double* yGradient;
    int width;
    int height;
    int pyramidDepth;
    Transformation transformation;
    bool isTarget;

    void init();

    void antiSymmetricFirMirrorOffBounds1D (
        double* h,
        int hlen,
        double* c,
        int clen,
        double* s,
        int slen
    );

    void basicToCardinal2D (
            double* basic,
            double* cardinal,
            int width,
            int height,
            int degree
    );
    void buildCoefficientPyramid ();
    void buildImageAndGradientPyramid ();
    void buildImagePyramid ();
    void cardinalToDual2D (
            double* cardinal,
            double* dual,
            int width,
            int height,
            int degree
    );

    /*------------------------------------------------------------------*/
    void coefficientToGradient1D (
            double* c,
            int clen
    );

    /*------------------------------------------------------------------*/
    void coefficientToSamples1D (
            double* c,
            int clen
    ) ;

    /*------------------------------------------------------------------*/
    void coefficientToXYGradient2D (
            double* basic,
            double* xGradient,
            double* yGradient,
            int width,
            int height
    );

    /*------------------------------------------------------------------*/
    void dualToCardinal2D (
            double* dual,
            double* cardinal,
            int width,
            int height,
            int degree
    );

    /*------------------------------------------------------------------*/
    void extractColumn (
            double* array,
            int width,
            int x,
            double* column,
            int height
    );

    /*------------------------------------------------------------------*/
    void extractRow (
            double* array,
            int y,
            double* row,
            int width
    );

    /*------------------------------------------------------------------*/
    double* getBasicFromCardinal2D (
    );

    /*------------------------------------------------------------------*/
    double* getBasicFromCardinal2D (
            double* cardinal,
            int width,
            int height,
            int degree
    );
    /*------------------------------------------------------------------*/
    double* getHalfDual2D (
            double* fullDual,
            int fullWidth,
            int fullHeight
    );

    /*------------------------------------------------------------------*/
    double getInitialAntiCausalCoefficientMirrorOffBounds (
            double* c,
            int clen,
            double z,
            double tolerance
    );

    /*------------------------------------------------------------------*/
    double getInitialCausalCoefficientMirrorOffBounds (
            double* c,
            int clen,
            double z,
            double tolerance
    );

    void imageToXYGradient2D (
    );

    /*------------------------------------------------------------------*/
    void putColumn (
            double* array,
            int width,
            int x,
            double* column,
            int columnlen
    );
    /*------------------------------------------------------------------*/
    void putRow (
        double* array,
        int y,
        double* row,
        int rowlen
    );
    /*------------------------------------------------------------------*/
    void reduceDual1D (
        double* c,
        int clen,
        double* s,
        int slen
    );
    
    void TurboRegImage::samplesToInterpolationCoefficient1D (
            double* c,
            int clen,
            int degree,
            double tolerance
    );

    void symmetricFirMirrorOffBounds1D (
        double* h,
        int hlen,
        double* c,
        int clen,
        double* s,
        int slen
    );

}
