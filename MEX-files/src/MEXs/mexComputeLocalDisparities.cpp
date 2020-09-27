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
    #define segL_IN prhs[0]
    #define numSegL_IN prhs [1]
    #define segR_IN prhs[2]
    #define numSegR_IN prhs [3]
    #define image1_IN prhs[4]
    #define image2_IN prhs[5]
    #define baseline_IN prhs[6]
    #define input_IN prhs[7]

    /* Macros for the output arguments */
    #define dispL_OUT plhs[0]
    #define scoresL_OUT plhs[1]
    #define dispR_OUT plhs[2]
    #define scoresR_OUT plhs[3]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 8 || nrhs > 8)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 4)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(segL_IN);
    int imgCols = mxGetN(segL_IN);

    int* segL = (int*)mxGetData(segL_IN);
    int numSegL = (int)mxGetScalar(numSegL_IN);
    int* segR = (int*)mxGetData(segR_IN);
    int numSegR = (int)mxGetScalar(numSegR_IN);

    unsigned char* image1 = (unsigned char*)mxGetData(image1_IN);
    unsigned char* image2 = (unsigned char*)mxGetData(image2_IN);

    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(image1_IN) > 2)
        imgChannels = 3;

    int baseline = (int)ceil(mxGetScalar(baseline_IN));

    // Input variables
    int maxDisparity = (int)(imgCols*mxGetScalar(mxGetProperty(input_IN, 0, "maxDisparity")));
    double dispLambda = mxGetScalar(mxGetProperty(input_IN, 0, "dispLambda"));
    // - Also for type of matching cost
    mxArray* cost = (mxGetProperty(input_IN, 0, "typeMatchCost"));
    int buflen = (mxGetM(cost) * mxGetN(cost)) + 1;
    char* mCost = (char *)mxCalloc(buflen, sizeof(char));
    mxGetString(cost, mCost, buflen);

    /* Create output data */
    const mwSize dims[2] = {imgRows, imgCols};
    dispL_OUT = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int* dispL = (int*)mxGetData(dispL_OUT);
    dispR_OUT = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int* dispR = (int*)mxGetData(dispR_OUT);

    const mwSize dims3[3] = {imgRows, imgCols, maxDisparity};
    scoresL_OUT = mxCreateNumericArray(3, dims3, mxDOUBLE_CLASS, mxREAL);
    double* scoresL = mxGetPr(scoresL_OUT);
    scoresR_OUT = mxCreateNumericArray(3, dims3, mxDOUBLE_CLASS, mxREAL);
    double* scoresR = mxGetPr(scoresR_OUT);

    /* Call method: computeQuickDisparities */
    computeLocalDisparities(dispL, scoresL, dispR, scoresR, imgRows, imgCols, segL, numSegL, segR, numSegR, image1, image2, imgChannels, maxDisparity, dispLambda, baseline, mCost);

    return;
}
