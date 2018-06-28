#include <iostream>
#include <vector>
#include <fstream>
#include <cstdint>
#include <iterator>

#include "TurboReg.h"
#include "TurboRegImage.h"
#include "TurboRegMask.h"
#include "TurboRegPointHandler.h"
#include "TurboRegTransform.h"
#include "matrix.h"

using namespace std;

std::vector<uint16_t> readFile(const char* filename)
{
   std::ifstream is;
   std::vector<uint16_t> rawfilebuffer;

   is.open(filename, std::ios::binary);
   is.seekg(0, std::ios::end);
   size_t filesize=is.tellg();
   is.seekg(0, std::ios::beg);

   rawfilebuffer.resize(filesize/sizeof(uint16_t));

   is.read((char *)rawfilebuffer.data(), filesize);

   return rawfilebuffer;

}

void writeFile(const char* filename, std::vector<uint16_t> &data)
{
   ofstream fout(filename, ios::out | ios::binary);
   fout.write((char*)&data[0], data.size() * sizeof(data[0]));
   fout.close();
}


int main(int argc, const char **argv)
{
    //int width = 128;
    //int height = 256;

	//swapped
	int width = 256;
	int height = 128;

    std::vector<uint16_t> int_imgdata_ref = readFile("/media/storage/eric/data/180614_Rut1_Papain/Use These Files/pygreg_test_ref.bin");
    std::vector<uint16_t> int_imgdata_mov = readFile("/media/storage/eric/data/180614_Rut1_Papain/Use These Files/pygreg_test_mov.bin");

    std::vector<double> imgdata_ref(int_imgdata_ref.begin(), int_imgdata_ref.end());
    std::vector<double> imgdata_mov(int_imgdata_mov.begin(), int_imgdata_mov.end());

    double *pDataRef= &imgdata_ref[0];
    double *pDataMov= &imgdata_mov[0];

    Transformation transformation = RIGID_BODY;

    TurboRegImage refImg(pDataRef, width, height, transformation, true);
    TurboRegImage movImg(pDataMov, width, height, transformation, false);

    TurboRegPointHandler refPH(refImg, transformation);
    TurboRegPointHandler movPH(movImg, transformation);

    TurboRegMask refMsk(refImg);
    TurboRegMask movMsk(movImg);

    refMsk.clearMask();
    movMsk.clearMask();

    int pyramidDepth = getPyramidDepth(
    		movImg.getWidth(), movImg.getHeight(),
			refImg.getWidth(), refImg.getHeight()
			);
    refImg.setPyramidDepth(pyramidDepth);
    refMsk.setPyramidDepth(pyramidDepth);
    movImg.setPyramidDepth(pyramidDepth);
    movMsk.setPyramidDepth(pyramidDepth);
    
    refImg.init();
    refMsk.init();
    movImg.init();
    movMsk.init();

  
    TurboRegTransform tform(&movImg, &movMsk, &movPH, &refImg, &refMsk, &refPH, transformation, false);
    
    tform.doRegistration();

    std::vector<double> imgout = tform.doFinalTransform(width, height);

    TurboRegPointHandler refPH2(refPH.getPoints());
    TurboRegPointHandler movPH2(movPH.getPoints());

    std::vector<double> imgout2 = tform.doFinalTransform (&movImg, &movPH2, &refImg, &refPH2, RIGID_BODY, false);

    matrix<double> tm = tform.getTransformationMatrix();

    //std::vector<double> imgout2 = tform.doFinalTransform(movImg, tm);

    tform.printPoints();

    std::vector<uint16_t> int_imgout(imgout2.begin(), imgout2.end());

    writeFile("/media/storage/eric/data/180614_Rut1_Papain/Use These Files/pygreg_test_out.bin", int_imgout);

	return 0;
}
