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
    #define dispL_IN prhs[0]
    #define dispR_IN prhs[1]

    /* Macros for the output arguments */
    #define disparity_OUT plhs[0]
    #define consistency_OUT plhs[1]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 2 || nrhs > 2)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(dispL_IN);
    int imgCols = mxGetN(dispL_IN);

    int* dispL = (int*)mxGetData(dispL_IN);
    int* dispR = (int*)mxGetData(dispR_IN);

    /* Create output data */
    const mwSize dims[2] = {imgRows, imgCols};
    disparity_OUT = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int* disparity = (int*)mxGetData(disparity_OUT);

    const mwSize dimsBool[2] = {imgRows, imgCols};
    consistency_OUT = mxCreateNumericArray(2, dimsBool, mxINT8_CLASS, mxREAL);
    bool* consistency = (bool*)mxGetData(consistency_OUT);

    /* Call method: crossValidation */
    crossValidation(disparity, consistency, imgRows, imgCols, dispL, dispR);

    return;
}
