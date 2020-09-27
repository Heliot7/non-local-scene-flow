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
    #define disparity_IN prhs[0]

    /* Macros for the output arguments */
    #define disparity_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 1 || nrhs > 1)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(disparity_IN);
    int imgCols = mxGetN(disparity_IN);

    double* disparity = mxGetPr(disparity_IN);

    /* Create output data */
    disparity_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* newDisparity = mxGetPr(disparity_OUT);

    /* Call method: consistentFill */
    consistentFill(newDisparity, imgRows, imgCols, disparity);

    return;
}
