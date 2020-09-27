#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Macros for the input arguments */
    #define projMatrices_IN prhs[0]
    #define focals_IN prhs[1]
    #define principals_IN prhs[2]
    #define camPos_IN prhs[3]
    #define Z_IN prhs[4]
    #define Z_x_IN prhs[5]
    #define Z_y_IN prhs[6]
    #define u_IN prhs[7]
    #define u_x_IN prhs[8]
    #define u_y_IN prhs[9]
    #define v_IN prhs[10]
    #define v_x_IN prhs[11]
    #define v_y_IN prhs[12]
    #define w_IN prhs[13]
    #define w_x_IN prhs[14]
    #define w_y_IN prhs[15]
    #define t0WarpedImgs_x_IN prhs[16]
    #define t0WarpedImgs_y_IN prhs[17]
    #define t1WarpedImgs_x_IN prhs[18]
    #define t1WarpedImgs_y_IN prhs[19]

    /* Macros for the output arguments */
    #define t0PartialUnknowns_OUT plhs[0]
    #define t1PartialUnknowns_OUT plhs[1]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 20 || nrhs > 20)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int numCameras = mxGetM(t0WarpedImgs_x_IN);
    int imgRows = mxGetM(Z_IN);
    int imgCols = mxGetN(Z_IN);
    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(t0WarpedImgs_x_IN) > 3)
        imgChannels = 3;

    double* projMatrices[numCameras];
    for(int CAMi = 0; CAMi < numCameras; ++CAMi)
        projMatrices[CAMi] = mxGetPr(mxGetCell(projMatrices_IN, CAMi));
    double* focals = mxGetPr(focals_IN);
    double* principals = mxGetPr(principals_IN);
    double* camPos = mxGetPr(camPos_IN);

    double* Z = mxGetPr(Z_IN);
    double* Z_x = mxGetPr(Z_x_IN);
    double* Z_y = mxGetPr(Z_y_IN);
    double* u = mxGetPr(u_IN);
    double* u_x = mxGetPr(u_x_IN);
    double* u_y = mxGetPr(u_y_IN);
    double* v = mxGetPr(v_IN);
    double* v_x = mxGetPr(v_x_IN);
    double* v_y = mxGetPr(v_y_IN);
    double* w = mxGetPr(w_IN);
    double* w_x = mxGetPr(w_x_IN);
    double* w_y = mxGetPr(w_y_IN);

    double* t0WarpedImgs_x = mxGetPr(t0WarpedImgs_x_IN);
    double* t0WarpedImgs_y = mxGetPr(t0WarpedImgs_y_IN);
    double* t1WarpedImgs_x = mxGetPr(t1WarpedImgs_x_IN);
    double* t1WarpedImgs_y = mxGetPr(t1WarpedImgs_y_IN);

    const mwSize dims[5] = {numCameras, imgRows, imgCols, imgChannels, 4}; // 4 unknowns: Z, u, v, w
    t0PartialUnknowns_OUT = mxCreateNumericArray(5, dims, mxDOUBLE_CLASS, mxREAL);
    double* t0PartialUnknowns = mxGetPr(t0PartialUnknowns_OUT);
    t1PartialUnknowns_OUT = mxCreateNumericArray(5, dims, mxDOUBLE_CLASS, mxREAL);
    double* t1PartialUnknowns = mxGetPr(t1PartialUnknowns_OUT);

    /* Call method: computePartialUnknowns */
    // Cases frame 't' & 't+1'
    computePartialUnknowns(t0PartialUnknowns, numCameras, imgRows, imgCols, imgChannels, projMatrices, focals, principals, camPos, Z, Z_x, Z_y, u, u_x, u_y, v, v_x, v_y, w, w_x, w_y, t0WarpedImgs_x, t0WarpedImgs_y, 0);
    computePartialUnknowns(t1PartialUnknowns, numCameras, imgRows, imgCols, imgChannels, projMatrices, focals, principals, camPos, Z, Z_x, Z_y, u, u_x, u_y, v, v_x, v_y, w, w_x, w_y, t1WarpedImgs_x, t1WarpedImgs_y, 1);
}

