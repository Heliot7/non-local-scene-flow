#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    /* Macros for the input arguments */
    #define t0Images_IN prhs[0]
    #define t1Images_IN prhs[1]
    #define t0ProjImgs_IN prhs[2]
    #define t1ProjImgs_IN prhs[3]
    #define t0occlusionMaps_IN prhs[4]
    #define t1occlusionMaps_IN prhs[5]

    /* Macros for the output arguments */
    #define t0WarpedImgs_OUT plhs[0]
    #define t1WarpedImgs_OUT plhs[1]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 6 || nrhs > 6)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(mxGetCell(t0Images_IN, 0));
    int imgCols = mxGetN(mxGetCell(t0Images_IN, 0));
    int numCameras = mxGetN(t0Images_IN);
    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(mxGetCell(t0Images_IN, 0)) > 2)
    {
        imgChannels = 3;
        imgCols /= 3;
    }

    double* t0ProjImgs = mxGetPr(t0ProjImgs_IN);
    double* t1ProjImgs = mxGetPr(t1ProjImgs_IN);
    bool* t0occlusionMaps = (bool*)mxGetData(t0occlusionMaps_IN);
    bool* t1occlusionMaps = (bool*)mxGetData(t1occlusionMaps_IN);

    unsigned char* t0Images[numCameras];
    unsigned char* t1Images[numCameras];
    for(int CAMi = 0; CAMi < numCameras; ++CAMi)
    {
        t0Images[CAMi] = (unsigned char*)mxGetData(mxGetCell(t0Images_IN, CAMi));
        t1Images[CAMi] = (unsigned char*)mxGetData(mxGetCell(t1Images_IN, CAMi));
    }

    const mwSize dims[4] = {numCameras, imgRows, imgCols, imgChannels};
    t0WarpedImgs_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double* t0WarpedImgs = mxGetPr(t0WarpedImgs_OUT);
    t1WarpedImgs_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double* t1WarpedImgs = mxGetPr(t1WarpedImgs_OUT);

    /* Call method: buildWarpedImage */
    buildWarpedImage(t0WarpedImgs, numCameras, imgRows, imgCols, imgChannels, t0Images, t0ProjImgs, t0occlusionMaps);
    buildWarpedImage(t1WarpedImgs, numCameras, imgRows, imgCols, imgChannels, t1Images, t1ProjImgs, t1occlusionMaps);

    return;
}
