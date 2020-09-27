#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Macros for the input arguments */
    #define dZ_IN prhs[0]
    #define du_IN prhs[1]
    #define dv_IN prhs[2]
    #define dw_IN prhs[3]
    #define t0WarpedImgs_IN prhs[4]
    #define t1WarpedImgs_IN prhs[5]
    #define t0I_unknown_IN prhs[6]
    #define t1I_unknown_IN prhs[7]
    #define input_IN prhs[8]

    /* Macros for the output arguments */
    #define psiData_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 9 || nrhs > 9)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int numCameras = mxGetM(t0WarpedImgs_IN);
    int imgRows = mxGetM(dZ_IN);
    int imgCols = mxGetN(dZ_IN);
    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(t0WarpedImgs_IN) > 3)
        imgChannels = 3;

    double* dZ = mxGetPr(dZ_IN);
    double* du = mxGetPr(du_IN);
    double* dv = mxGetPr(dv_IN);
    double* dw = mxGetPr(dw_IN);

    double* t0WarpedImgs = mxGetPr(t0WarpedImgs_IN);
    double* t1WarpedImgs = mxGetPr(t1WarpedImgs_IN);
    double* t0I_unknown = mxGetPr(t0I_unknown_IN);
    double* t1I_unknown = mxGetPr(t1I_unknown_IN);

    double epsilon = mxGetScalar(mxGetProperty(input_IN, 0, "epsData"));

    const mwSize dims[5] = {numCameras, imgRows, imgCols, imgChannels, 3}; // 3 types of connections (2 spatial and 1 temporal)
    psiData_OUT = mxCreateNumericArray(5, dims, mxDOUBLE_CLASS, mxREAL);
    double* psiData = mxGetPr(psiData_OUT);

    /* Call method: buildRobustnessData */
    buildRobustnessData(psiData, numCameras, imgRows, imgCols, imgChannels, dZ, du, dv, dw, t0WarpedImgs, t1WarpedImgs, t0I_unknown, t1I_unknown, epsilon);
}

