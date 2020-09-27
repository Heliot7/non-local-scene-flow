#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libImgProc/segmentation.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Macros for the input arguments */
    #define image_IN prhs[0]
    #define Z_IN prhs[1]
    #define u_IN prhs [2]
    #define v_IN prhs [3]
    #define w_IN prhs [4]
    #define input_IN prhs[5]

    /* Macros for the output arguments */
    #define segmentation_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 6 || nrhs > 6)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    mwSize imgRows = mxGetM(image_IN);
    mwSize imgCols = mxGetN(image_IN);
    unsigned char* image = (unsigned char*)mxGetData(image_IN);
    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(image_IN) > 2)
    {
        imgChannels = 3;
        imgCols /= 3;
    }

    double* Z = mxGetPr(Z_IN);
    double* u = mxGetPr(u_IN);
    double* v = mxGetPr(v_IN);
    double* w = mxGetPr(w_IN);

    // Input variables
    InputKMeans input;
    input.numK = mxGetScalar(mxGetProperty(input_IN, 0, "numK"));
    input.numIterKMeans = mxGetScalar(mxGetProperty(input_IN, 0, "numIterKMeans"));
    input.onConnectedComponents = mxGetScalar(mxGetProperty(input_IN, 0, "onConnectedComponents"));

    /* Create output data */    
    const mwSize dims[2] = {imgRows, imgCols};
    segmentation_OUT = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int* segmentation = (int*)mxGetData(segmentation_OUT);

    /* Call method: KMeans9D */
    KMeans9D(segmentation, imgRows, imgCols, image, imgChannels, Z, u, v, w, input);

    return;
}
