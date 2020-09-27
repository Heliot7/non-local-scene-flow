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
    #define segmentation_IN prhs[1]
    #define Z_IN prhs[2]
    #define u_IN prhs [3]
    #define v_IN prhs [4]
    #define w_IN prhs [5]
    #define mergedThreshold_IN prhs[6]

    /* Macros for the output arguments */
    #define segmentation_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 7 || nrhs > 7)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int* segmentation = (int*)mxGetData(segmentation_IN);
    mwSize imgRows = mxGetM(segmentation_IN);
    mwSize imgCols = mxGetN(segmentation_IN);
    unsigned char* image = (unsigned char*)mxGetData(image_IN);

    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(image_IN) > 2)
        imgChannels = 3;

    double* Z = mxGetPr(Z_IN);
    double* u = mxGetPr(u_IN);
    double* v = mxGetPr(v_IN);
    double* w = mxGetPr(w_IN);

    double mergedTh = mxGetScalar(mergedThreshold_IN);

    /* Create output data */
    const mwSize dims[2] = {imgRows, imgCols};
    segmentation_OUT = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int* segmentsOut = (int*)mxGetData(segmentation_OUT);

    int imgLabels = 0;
    for(int i = 0; i < imgRows*imgCols; ++i)
    {
        segmentsOut[i] = segmentation[i];
        if(segmentsOut[i] > imgLabels)
            imgLabels = segmentsOut[i];
    }

    for(int i = 0; i < imgRows*imgCols; ++i)
        segmentsOut[i]--;

    /* Call method: ConnectedSegments */
    mergeSegments(imgRows, imgCols, imgChannels, imgLabels, segmentsOut, image, Z, u, v, w, mergedTh);

    for(int i = 0; i < imgRows*imgCols; ++i)
        segmentsOut[i]++;

    return;
}
