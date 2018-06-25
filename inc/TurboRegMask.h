#include <stack>

//using namespace std;

class TurboRegMask
{ 

public:
    TurboRegMask (double* imp, int width, int height);
    
    void clearMask ();
    double* getMask ();

    std::stack <double*> getPyramid();

    void setPyramidDepth (int pyramidDepth);

    virtual ~TurboRegMask();

private:

    std::stack<double*> pyramid;

    double* mask;
    int width;
    int height;
    int pyramidDepth;

    void init ();
    void buildPyramid ();

    double* getHalfMask2D (
            double* fullMask,
            int fullWidth,
            int fullHeight
    );
        

} /* end class turboRegMask */
