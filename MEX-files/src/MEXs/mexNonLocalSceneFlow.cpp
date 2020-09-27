#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libImgProc/filters.hpp"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    /* Macros for the input arguments */
    #define matrix_IN prhs[0]
    #define windowSize_IN prhs[1]

    /* Macros for the output arguments */
    #define matrixMedian_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 2 || nrhs > 2)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    mwSize imgRows = mxGetM(matrix_IN);
    mwSize imgCols = mxGetN(matrix_IN);
    double* matrix = mxGetPr(matrix_IN);
    int windowSize = mxGetScalar(windowSize_IN);

    /* Create output data */
    matrixMedian_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* matrixMedian = mxGetPr(matrixMedian_OUT);

    /* Call method: KMeans9D */
    medianFilter(matrixMedian, imgRows, imgCols, matrix, windowSize);

    return;
}
