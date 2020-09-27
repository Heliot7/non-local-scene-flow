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
    #define t0ProjImgs_IN prhs[0]
    #define t1ProjImgs_IN prhs[1]
    #define t0occlusionMaps_IN prhs[2]
    #define t1occlusionMaps_IN prhs[3]

    /* Macros for the output arguments */
    #define occlusionBool_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 4 || nrhs > 4)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int numCameras = mxGetM(t0occlusionMaps_IN);
    int imgRows = mxGetDimensions(t0occlusionMaps_IN)[1];
    int imgCols = mxGetDimensions(t0occlusionMaps_IN)[2];

    double* t0ProjImgs = mxGetPr(t0ProjImgs_IN);
    double* t1ProjImgs = mxGetPr(t1ProjImgs_IN);
    bool* t0occlusionMaps = (bool*)mxGetData(t0occlusionMaps_IN);
    bool* t1occlusionMaps = (bool*)mxGetData(t1occlusionMaps_IN);

    const mwSize dims[4] = {numCameras, imgRows, imgCols, 3}; // 3 types of connections (2 spatial and 1 temporal)
    occlusionBool_OUT = mxCreateNumericArray(4, dims, mxUINT8_CLASS, mxREAL);
    bool* occBin = (bool*)mxGetData(occlusionBool_OUT);

    /* Call method: buildOcclusionBinary */
    buildOcclusionBinary(occBin, numCameras, imgRows, imgCols, t0ProjImgs, t1ProjImgs, t0occlusionMaps, t1occlusionMaps);

    return;
}
