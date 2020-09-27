#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libImgProc/disparity.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Macros for the input arguments */
    #define segmentation_IN prhs[0]
    #define numSegments_IN prhs [1]
    #define planes_IN prhs [2]
    #define disparity_IN prhs[3]
    #define consistency_IN prhs[4]
    #define input_IN prhs[5]
    #define LEVEL_IN prhs[6]

    /* Macros for the output arguments */
    #define similarities_OUT plhs[0]

    // Number of input/output arguments.
    if(nrhs < 7 || nrhs > 7)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(segmentation_IN);
    int imgCols = mxGetN(segmentation_IN);

    int* segmentation = (int*)mxGetData(segmentation_IN);
    int numSegments = (int)mxGetScalar(numSegments_IN);
    double* disparity = mxGetPr(disparity_IN);
    bool* consistency = (bool*)mxGetData(consistency_IN);
    double* planes = mxGetPr(planes_IN);

    // Input variables
    double pImpact = mxGetScalar(mxGetProperty(input_IN, 0, "pImpact"));
    int LEVEL = (int)mxGetScalar(LEVEL_IN);
    int numLevel = mxGetScalar(mxGetProperty(input_IN, 0, "endLevel")) - LEVEL;
    bool onImpact = (bool)mxGetScalar(mxGetProperty(input_IN, 0, "onImpact"));
    if(onImpact)
        pImpact *= exp(-numLevel);
    mexPrintf("pImpact %f at LEVEL %d\n", pImpact, LEVEL);

    /* Create output data */
    const mwSize dims[3] = {imgRows, imgCols, 4};
    similarities_OUT = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* similarities = mxGetPr(similarities_OUT);

    /* Call method: computeSegmentPlanes */
    planeSimilarities(similarities, imgRows, imgCols, segmentation, numSegments, planes, disparity, consistency, pImpact);

    return;
}
