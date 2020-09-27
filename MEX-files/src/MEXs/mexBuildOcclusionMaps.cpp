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
    #define t0ProjImgs_IN prhs[0]
    #define t1ProjImgs_IN prhs[1]
    #define X_IN prhs[2]
    #define Y_IN prhs[3]
    #define Z_IN prhs[4]
    #define u_IN prhs[5]
    #define v_IN prhs[6]
    #define w_IN prhs[7]
    #define COP_IN prhs[8]

    /* Macros for the output arguments */
    #define t0occlusionMaps_OUT plhs[0]
    #define t1occlusionMaps_OUT plhs[1]

    /* Check correctness of input/output arguments */
    // Number of input/output arguments.
    if(nrhs < 9 || nrhs > 9)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(X_IN);
    int imgCols = mxGetN(X_IN);
    int numCameras = mxGetM(t0ProjImgs_IN);

    double* t0ProjImgs = mxGetPr(t0ProjImgs_IN);
    double* t1ProjImgs = mxGetPr(t1ProjImgs_IN);
    double* X = mxGetPr(X_IN);
    double* Y = mxGetPr(Y_IN);
    double* Z = mxGetPr(Z_IN);
    double* u = mxGetPr(u_IN);
    double* v = mxGetPr(v_IN);
    double* w = mxGetPr(w_IN);

    double* COP[numCameras];
    for(int CAMi = 0; CAMi < numCameras; ++CAMi)
        COP[CAMi] = mxGetPr(mxGetCell(COP_IN, CAMi));

    const mwSize dims[3] = {numCameras, imgRows, imgCols};
    t0occlusionMaps_OUT = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    bool* t0occlusionMaps = (bool*)mxGetData(t0occlusionMaps_OUT);
    t1occlusionMaps_OUT = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
    bool* t1occlusionMaps = (bool*)mxGetData(t1occlusionMaps_OUT);

    /* Call method: buildOcclusionMaps */
    buildOcclusionMaps(t0occlusionMaps, numCameras, imgRows, imgCols, t0ProjImgs, X, Y, Z, u, v, w, COP);
    buildOcclusionMaps(t1occlusionMaps, numCameras, imgRows, imgCols, t1ProjImgs, X, Y, Z, u, v, w, COP);

    return;
}
