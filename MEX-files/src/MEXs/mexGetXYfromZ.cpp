#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    /* Macros for the input arguments */
    #define Z_IN prhs[0]
    #define focals_IN prhs[1]
    #define principals_IN prhs[2]
    #define COP_IN prhs[3]

    /* Macros for the output arguments */
    #define X_OUT plhs[0]
    #define Y_OUT plhs[1]

    /* Check correctness of input/output arguments */
    // Number of input/output arguments.
    if(nrhs < 4 || nrhs > 4)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(Z_IN);
    int imgCols = mxGetN(Z_IN);
    double* Z = mxGetPr(Z_IN);

    double* focals = mxGetPr(focals_IN);
    double* principals = mxGetPr(principals_IN);
    double* COP = mxGetPr(COP_IN);

    /* Call method: getXYfromZ */
    X_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* X = mxGetPr(X_OUT);
    Y_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* Y = mxGetPr(Y_OUT);
    getXYfromZ(imgRows, imgCols, X, Y, Z, focals, principals, COP);

    return;
}
