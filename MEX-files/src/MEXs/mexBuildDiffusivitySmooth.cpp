#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>

#include "mex.h"
#include "engine.h"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Macros for the input arguments */
    #define Z_x_IN prhs[0]
    #define Z_y_IN prhs[1]
    #define u_x_IN prhs[2]
    #define u_y_IN prhs[3]
    #define v_x_IN prhs[4]
    #define v_y_IN prhs[5]
    #define w_x_IN prhs[6]
    #define w_y_IN prhs[7]
    #define dZ_x_IN prhs[8]
    #define dZ_y_IN prhs[9]
    #define du_x_IN prhs[10]
    #define du_y_IN prhs[11]
    #define dv_x_IN prhs[12]
    #define dv_y_IN prhs[13]
    #define dw_x_IN prhs[14]
    #define dw_y_IN prhs[15]
    #define input_IN prhs[16]

    /* Macros for the output arguments */
    #define diffSmooth_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 17 || nrhs > 17)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    /* Get input data */
    int imgRows = mxGetM(Z_x_IN);
    int imgCols = mxGetN(Z_x_IN);
    double* Z_x = mxGetPr(Z_x_IN);
    double* Z_y = mxGetPr(Z_y_IN);
    double* u_x = mxGetPr(u_x_IN);
    double* u_y = mxGetPr(u_y_IN);
    double* v_x = mxGetPr(v_x_IN);
    double* v_y = mxGetPr(v_y_IN);
    double* w_x = mxGetPr(w_x_IN);
    double* w_y = mxGetPr(w_y_IN);
    double* dZ_x = mxGetPr(dZ_x_IN);
    double* dZ_y = mxGetPr(dZ_y_IN);
    double* du_x = mxGetPr(du_x_IN);
    double* du_y = mxGetPr(du_y_IN);
    double* dv_x = mxGetPr(dv_x_IN);
    double* dv_y = mxGetPr(dv_y_IN);
    double* dw_x = mxGetPr(dw_x_IN);
    double* dw_y = mxGetPr(dw_y_IN);

    double epsilon = mxGetScalar(mxGetProperty(input_IN, 0, "epsSmooth"));
    double muZ = mxGetScalar(mxGetProperty(input_IN, 0, "muZ"));
    double muUVW = mxGetScalar(mxGetProperty(input_IN, 0, "muUVW"));


    const mwSize dims[3] = {imgRows, imgCols, 2}; // 2 kinds of diffusivity: space (Z) and motion (uvw)
    diffSmooth_OUT = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* diffSmooth = mxGetPr(diffSmooth_OUT);

    /* Call method: buildDiffusivitySmooth */
    buildDiffusivitySmooth(diffSmooth, imgRows, imgCols, Z_x, Z_y, u_x, u_y, v_x, v_y, w_x, w_y, dZ_x, dZ_y, du_x, du_y, dv_x, dv_y, dw_x, dw_y, epsilon, muZ, muUVW);
}

