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
    #define disparity_IN prhs[2]
    #define scores_IN prhs[3]
    #define consistency_IN prhs[4]
    #define input_IN prhs[5]

    /* Macros for the output arguments */
    #define disparity_OUT plhs[0]
    #define planes_OUT plhs[1]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 6 || nrhs > 6)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(segmentation_IN);
    int imgCols = mxGetN(segmentation_IN);

    int* segmentation = (int*)mxGetData(segmentation_IN);
    int numSegments = (int)mxGetScalar(numSegments_IN);
    int* disparity = (int*)mxGetData(disparity_IN);
    double* scores = mxGetPr(scores_IN);
    bool* consistency = (bool*)mxGetData(consistency_IN);

    // Input variables
    int maxDisparity = (int)(imgCols*mxGetScalar(mxGetProperty(input_IN, 0, "maxDisparity")));
    int maxIterPlane = mxGetScalar(mxGetProperty(input_IN, 0, "maxIterPlane"));

    /* Create output data */
    disparity_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* newDisparity = mxGetPr(disparity_OUT);
    planes_OUT = mxCreateDoubleMatrix(numSegments, 3, mxREAL);
    double* planes = mxGetPr(planes_OUT);

    /* Call method: computeSegmentPlanes */
    computeSegmentPlanes(newDisparity, planes, imgRows, imgCols, segmentation, numSegments, disparity, scores, consistency, maxDisparity, maxIterPlane);

    return;
}
