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
    #define a_IN prhs[0]
    #define b_IN prhs[1]
    #define psiSmooth_IN prhs[2]
    #define segmentPenalties_IN prhs [3]
    #define LEVEL_IN prhs[4]
    #define input_IN prhs[5]

    /* Macros for the output arguments */
    #define NEWdZ_OUT plhs[0]
    #define NEWdu_OUT plhs[1]
    #define NEWdv_OUT plhs[2]
    #define NEWdw_OUT plhs[3]
    #define evolResidual_OUT plhs[4]

    /* Check correctness of input/output arguments */
    // Number of input/output arguments.
    if(nrhs < 6 || nrhs > 6)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 5)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetDimensions(a_IN)[0];
    int imgCols = mxGetDimensions(a_IN)[1];

    double* a = mxGetPr(a_IN);
    double* b = mxGetPr(b_IN);
    double* psiSmooth = mxGetPr(psiSmooth_IN);
    double* segmentPenalties = mxGetPr(segmentPenalties_IN);

    int LEVEL = (int)mxGetScalar(LEVEL_IN);
    int maxIter = mxGetScalar(mxGetProperty(input_IN, 0, "maxIter"));
    double errorSolver = mxGetScalar(mxGetProperty(input_IN, 0, "errorSolver"));
    double wSOR = mxGetScalar(mxGetProperty(input_IN, 0, "weightSOR"));
    int numLevels = (int)mxGetScalar(mxGetProperty(input_IN, 0, "numLevels"));
    double factor = mxGetScalar(mxGetProperty(input_IN, 0, "factor"));
    Weights weights;
    weights.alphaUV = mxGetScalar(mxGetProperty(input_IN, 0, "alphaUV"));
    weights.alphaW = mxGetScalar(mxGetProperty(input_IN, 0, "alphaW"));
    weights.alphaZ = mxGetScalar(mxGetProperty(input_IN, 0, "alphaZ"));
    weights.muUVW = mxGetScalar(mxGetProperty(input_IN, 0, "muUVW"));
    weights.muZ = mxGetScalar(mxGetProperty(input_IN, 0, "muZ"));

    // Get string with the used iterative solver.
    mxArray* solverMatlab = (mxGetProperty(input_IN, 0, "itSolver"));
    int buflen = (mxGetM(solverMatlab) * mxGetN(solverMatlab)) + 1;
    char* itSolver = (char *)mxCalloc(buflen, sizeof(char));
    mxGetString(solverMatlab, itSolver, buflen);
    // And preconditioner.
    mxArray* preconditioner = (mxGetProperty(input_IN, 0, "preconditioner"));
    buflen = (mxGetM(preconditioner) * mxGetN(preconditioner)) + 1;
    char* itPrecon = (char *)mxCalloc(buflen, sizeof(char));
    mxGetString(preconditioner, itPrecon, buflen);

    NEWdZ_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* dZ = mxGetPr(NEWdZ_OUT);
    NEWdu_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* du = mxGetPr(NEWdu_OUT);
    NEWdv_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* dv = mxGetPr(NEWdv_OUT);
    NEWdw_OUT = mxCreateDoubleMatrix(imgRows, imgCols, mxREAL);
    double* dw = mxGetPr(NEWdw_OUT);

    evolResidual_OUT = mxCreateDoubleMatrix(maxIter, 1, mxREAL);
    double* evolResidual = mxGetPr(evolResidual_OUT);
    for(int i = 0; i < maxIter; ++i)
        evolResidual[i] = -1.0;

    /* Call method: iterativeSolver */
    iterativeSolver(dZ, du, dv, dw, evolResidual, itSolver, itPrecon, imgRows, imgCols, a, b, psiSmooth, wSOR, maxIter, errorSolver, LEVEL, numLevels, factor, weights, segmentPenalties);
}


