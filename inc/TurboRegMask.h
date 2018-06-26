#include <stack>
#include "matrix.h"

#include "TurboRegImage.h"

#ifndef TURBOREGMASK_H_
#define TURBOREGMASK_H_

class MaskStackItem {
public:
    std::vector<double> halfMask;

    MaskStackItem(int size) {
        halfMask.resize(size);
    }
};

class TurboRegMask
{ 

public:
    TurboRegMask (TurboRegImage &img);
    TurboRegMask (matrix<double> &imp, int width, int height);
    
    void clearMask ();
    std::vector<double> &getMask();

    std::stack<MaskStackItem> getPyramid();

    void setPyramidDepth (int pyramidDepth);
    void init();

private:

    std::stack<MaskStackItem> pyramid;

    std::vector<double> mask;
    int width;
    int height;
    int pyramidDepth;

    void buildPyramid ();

    std::vector<double> getHalfMask2D (
            std::vector<double> &fullMask,
            int fullWidth,
            int fullHeight,
            std::vector<double> &halfMask
    );
        

}; /* end class turboRegMask */

#endif
