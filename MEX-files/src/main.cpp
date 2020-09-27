#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <unistd.h>
#include <typeinfo>

#include "libIO/file.hpp"
#include "libIO/lodepng.hpp"
#include "libIO/pgm.hpp"
#include "libIO/viewMatlab.hpp"

#include "libMath/metric.hpp"
#include "libMath/matrix.hpp"

#include "libSceneFlow/nonLocalSceneFlow.hpp"

#include "engine.h"

using namespace std;

// .. -DCMAKE_BUILD_TYPE=Debug <- CMake Compile arguments

int main()
{   
    /*
    // Check Matrix function
    Matrix<int> Mat(3,3);

    // Init input matrix
    Mat(0,0) =  1; Mat(0,1) = 3; Mat(0,2) =  5;
    Mat(1,0) = -1; Mat(1,1) = 2; Mat(1,2) =  0;
    Mat(2,0) =  4; Mat(2,1) = 2; Mat(2,2) = -3;

    cout << "Determinant 3x3: " << Mat.det3() << endl;

    // Solve inverse
    Matrix<double> invMat;
    Matrix<double>::inverse3(invMat, Mat);

    // Check solution
    Mat.print("Mat");
    invMat.print("Inverse Mat");

    return 0;
    */

    // -> Declare input parameters
    Input input;
    // [GLOBAL PARAMETERS]
    input.onSequence = false;
    input.dataset = "Teddy";
    input.imgExt = "png";
    input.numCameras = 2;
    input.onStore = false;
    input.onDisplay = false;
    input.onResidual = false;
    // [SCENE FLOW PARAMETERS]
    input.numLevels = 10;
    input.initLevel = 5;
    input.endLevel = 5;
    input.maxOut = 5;
    input.maxIn = 5;
    input.epsOut = -0.001;
    input.epsIn = -0.001;
    input.Zinit = 500;
    input.onZinit = false;
    input.border = 1;
    input.onSmoothing = true;
    input.alphaZ = 2.0;
    input.alphaUV = 30.0;
    input.alphaW = 30.0;
    input.mu = 1.0;
    input.muZ = 1.0;
    input.muUVW = 1.0;
    input.itSolver = "SOR";
    input.preconditioner = "NO";
    input.ILUp = 0;
    input.maxIter = 20;
    input.weightSOR = 1.8;
    input.errorSolver = 1e-2;
    input.epsSmooth = 0.0001;
    input.epsData = 0.0001;
    // [SEGMENTATION BASED]
    input.onSegmentationBased = false;
    input.proportionKmax = 0.1;
    input.distanceKmax = 0.1;
    input.numBins = 100;
    input.numIterKMeans = 10;
    input.maxTry = 5;
    input.maxDisparity = 0.35;
    input.maxIterPlane = 5;
    input.pImpact = 2.5;
    input.aggregatedWindowSize = 0.05;
    input.dispLambda = 0.01;

    // -> Path of image files
    string path;
    if(input.onSequence)
        path = "/work/panareda/data/";
    else
        path = "../../../data/";

    path.append(input.dataset.c_str());
    path.append("/");
    path.append(input.dataset.c_str());
    path.append("_t0_0.");
    path.append(input.imgExt.c_str());
    cout << path.c_str() << endl;
    if( fexists( path.c_str() ) )
        cout << "Input files exists!" << endl;

    /*
    // -> Loading all camera parameters
    CameraParameters infoCameras[input.numCameras];

    // -> Loop over all sequence of images
    int currentFrame = 0;
    bool goSequence = true;

    while(goSequence)
    {
        double* img0[input.numCameras];
        string imagePath;
        if( fexists( path.c_str() ) )
            cout << "Input files exists!" << endl;
        else
            cout << "Error: Input files do not exist!" << endl;


    }
    */

    unsigned width, height;
    int nChannels;
    unsigned char* image;

    //NonLocalSceneFlow* sceneFlow = new NonLocalSceneFlow(input);
    //sceneFlow->run();
    //return 0;

    // PNG files
    if(input.imgExt == "png")
    {
        cout << "Loading PNG file..." << endl;
        std::vector<unsigned char> buffer;
        std::vector<unsigned char> imageVec;
        lodepng::State state;
        lodepng::load_file(buffer, path.c_str());
        lodepng::decode(imageVec, width, height, state, buffer);
        LodePNGColorMode& color = state.info_png.color;
        nChannels = lodepng_get_channels(&color);
        cout << "num channels " << nChannels << endl;
        image = new unsigned char[width*height*nChannels];
        for(int i = 0; i < width; ++i)              // COL
            for(int j = 0; j < height; ++j)         // ROW
                for(int c = 0; c < nChannels; ++c)  // CHANNEL
                    image[height*i + j + c*height*width] = imageVec.at(4*width*j + 4*i + c+1); // 4 Channels in vector RGBA (TODO REVISE!!)
    }

    // PGM files
    PGMData *data;
    if(input.imgExt == "pgm")
    {
        cout << "Loading PGM file..." << endl;
        readPGM(path.c_str(), data);
        width = data->col;
        height = data->row;
        nChannels = 1;
        image = new unsigned char[width*height];
        for(int i = 0; i < width; ++i)
            for(int j = 0; j < height; ++j)
                image[height*i + j] = (unsigned char)data->matrix[j][i];
    }


    // NEEED TO CONVERTIT INTO
    cout << "File read!" << endl;
    showImage(image, height, width, nChannels);
    cout << "END" << endl;

    return 0;
}
