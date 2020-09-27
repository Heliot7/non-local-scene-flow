#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libImgProc/differentiation.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    /* Macros for the input arguments */
    #define matrix_IN prhs[0]

    /* Macros for the output arguments */
    #define dX_OUT plhs[0]
    #define dY_OUT plhs[1]

    /* Check correctness of input/output arguments */
    // Number of input/output arguments.
    if(nrhs < 1 || nrhs > 1)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(matrix_IN);
    int imgCols = mxGetN(matrix_IN);
    double* matrix = mxGetPr(matrix_IN);

    /* Call method: getXYfromZ */
    dX_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* dX = mxGetPr(dX_OUT);
    dY_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* dY = mxGetPr(dY_OUT);
    gradientDouble(imgRows, imgCols, matrix, dX, dY);

    return;
}
