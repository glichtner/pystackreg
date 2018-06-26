#include "TurboReg.h"

int getPyramidDepth (
        int sw,
        int sh,
        int tw,
        int th
) {
    int pyramidDepth = 1;
    while (((2 * PYRAMID_MIN_SIZE) <= sw)
            && ((2 * PYRAMID_MIN_SIZE) <= sh)
            && ((2 * PYRAMID_MIN_SIZE) <= tw)
            && ((2 * PYRAMID_MIN_SIZE) <= th)) {
        sw /= 2;
        sh /= 2;
        tw /= 2;
        th /= 2;
        pyramidDepth++;
    }
    return(pyramidDepth);
} /* end getPyramidDepth */
