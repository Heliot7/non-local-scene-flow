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
    #define Z_IN prhs[0]
    #define u_IN prhs[1]
    #define v_IN prhs[2]
    #define w_IN prhs[3]
    #define t0WarpedImgs_IN prhs[4]
    #define t1WarpedImgs_IN prhs[5]
    #define t0I_unknown_IN prhs[6]
    #define t1I_unknown_IN prhs[7]
    #define psiData_IN prhs[8]
    #define psiSmooth_IN prhs[9]
    #define segmentPenalties_IN prhs[10]
    #define occlusionBinary_IN prhs[11]
    #define LEVEL_IN prhs[12]
    #define input_IN prhs[13]

    /* Macros for the output arguments */
    #define matrixA_OUT plhs[0]
    #define matrixB_OUT plhs[1]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 14 || nrhs > 14)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int numCameras = mxGetM(t0WarpedImgs_IN);
    int imgRows = mxGetM(Z_IN);
    int imgCols = mxGetN(Z_IN);
    // Greyscale scheme
    int imgChannels = 1;
    // RGB scheme
    if(mxGetNumberOfDimensions(t0WarpedImgs_IN) > 3)
        imgChannels = 3;

    double* Z = mxGetPr(Z_IN);
    double* u = mxGetPr(u_IN);
    double* v = mxGetPr(v_IN);
    double* w = mxGetPr(w_IN);

    double* t0WarpedImgs = mxGetPr(t0WarpedImgs_IN);
    double* t1WarpedImgs = mxGetPr(t1WarpedImgs_IN);
    double* t0I_unknown = mxGetPr(t0I_unknown_IN);
    double* t1I_unknown = mxGetPr(t1I_unknown_IN);
    double* psiData = mxGetPr(psiData_IN);
    double* psiSmooth = mxGetPr(psiSmooth_IN);
    double* segmentPenalties = mxGetPr(segmentPenalties_IN);
    bool* occlusionBinary = (bool*)mxGetData(occlusionBinary_IN);

    int LEVEL = mxGetScalar(LEVEL_IN);
    int numLevels = (int)mxGetScalar(mxGetProperty(input_IN, 0, "numLevels"));
    double factor = mxGetScalar(mxGetProperty(input_IN, 0, "factor"));
    Weights weights;
    weights.alphaUV = mxGetScalar(mxGetProperty(input_IN, 0, "alphaUV"));
    weights.alphaW = mxGetScalar(mxGetProperty(input_IN, 0, "alphaW"));
    weights.alphaZ = mxGetScalar(mxGetProperty(input_IN, 0, "alphaZ"));
    weights.muUVW = mxGetScalar(mxGetProperty(input_IN, 0, "muUVW"));
    weights.muZ = mxGetScalar(mxGetProperty(input_IN, 0, "muZ"));

    const mwSize dimsA[4] = {imgRows, imgCols, 4, 4}; // 4x4 matrix per pixel
    matrixA_OUT = mxCreateNumericArray(4, dimsA, mxDOUBLE_CLASS, mxREAL);
    double* mA = mxGetPr(matrixA_OUT);
    const mwSize dimsB[3] = {imgRows, imgCols, 4}; // 4x4 matrix per pixel
    matrixB_OUT = mxCreateNumericArray(3, dimsB, mxDOUBLE_CLASS, mxREAL);
    double* mB = mxGetPr(matrixB_OUT);

    /* Call method: buildMatricesAb */
    buildMatricesAb(mA, mB, numCameras, imgRows, imgCols, imgChannels, Z, u, v, w, t0WarpedImgs, t1WarpedImgs,
                                      t0I_unknown, t1I_unknown, psiData, psiSmooth, occlusionBinary, LEVEL, numLevels, factor, weights, segmentPenalties);
}

