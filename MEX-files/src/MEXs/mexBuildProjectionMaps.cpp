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
    #define projMatrices_IN prhs[0]
    #define X_IN prhs[1]
    #define Y_IN prhs[2]
    #define Z_IN prhs[3]
    #define u_IN prhs[4]
    #define v_IN prhs[5]
    #define w_IN prhs[6]

    /* Macros for the output arguments */
    #define t0ProjImgs_OUT plhs[0]
    #define t1ProjImgs_OUT plhs[1]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 7 || nrhs > 7)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(X_IN);
    int imgCols = mxGetN(X_IN);
    int numCameras = mxGetN(projMatrices_IN);

    double* projMatrices[numCameras];
    for(int CAMi = 0; CAMi < numCameras; ++CAMi)
        projMatrices[CAMi] = mxGetPr(mxGetCell(projMatrices_IN, CAMi));

    double* X = mxGetPr(X_IN);
    double* Y = mxGetPr(Y_IN);
    double* Z = mxGetPr(Z_IN);
    double* u = mxGetPr(u_IN);
    double* v = mxGetPr(v_IN);
    double* w = mxGetPr(w_IN);

    const mwSize dims[4] = {numCameras, imgRows, imgCols, 2};
    t0ProjImgs_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double* t0ProjMap = mxGetPr(t0ProjImgs_OUT);
    t1ProjImgs_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    double* t1ProjMap = mxGetPr(t1ProjImgs_OUT);

    /* Call method: buildProjectionMaps */
    buildProjectionMaps(t0ProjMap, imgRows, imgCols, numCameras, projMatrices, X, Y, Z);
    buildProjectionMaps(t1ProjMap, imgRows, imgCols, numCameras, projMatrices, X, Y, Z, u, v, w);

    return;
}
