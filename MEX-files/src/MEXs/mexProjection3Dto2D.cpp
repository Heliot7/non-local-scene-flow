#include <iostream>
#include <stdio.h>
#include <math.h>

#include "mex.h"
#include "libSceneFlow/nonLocalSceneFlow.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    /* Macros for the input arguments */
    #define projMatrix_IN prhs[0]
    #define P_IN prhs[1]

    /* Macros for the output arguments */
    #define p_OUT plhs[0]

    /* Check correctness of input/output arguments */

    // Number of input/output arguments.
    if(nrhs < 2 || nrhs > 2)
        mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

    // Check type of input arguments
    if(!mxIsDouble(projMatrix_IN))
        mexErrMsgTxt("projMatrix must be a real 2D double array.");
    if(!mxIsDouble(P_IN))
        mexErrMsgTxt("P Matrix must be a real 1D double array.");

    /* Get input data */
    double* projMatrix = mxGetPr(projMatrix_IN);
    if ( mxGetM(projMatrix_IN)!= 3 || mxGetN(projMatrix_IN)!= 4 )
    {
        mexPrintf("matrix projMatrix is %dx%d\n", mxGetM(projMatrix_IN), mxGetN(projMatrix_IN));
        mexErrMsgTxt("projMatrix dimension is wrong (it should be 3x4)");
    }
    double* P = mxGetPr(P_IN);
    if ( mxGetM(P_IN)!= 3 || mxGetN(P_IN)!= 1 )
    {
        mexPrintf("matrix P is %dx%d\n", mxGetM(P_IN), mxGetN(P_IN));
        mexErrMsgTxt("P dimension is wrong (it should be 3x1)");
    }

    /* Call method: projection3Dto2D */
    P2D p = projection3Dto2D(projMatrix, P);

    /* Create output data */
    p_OUT = mxCreateNumericMatrix(2, 1, mxDOUBLE_CLASS, mxREAL);
    double* aPr = mxGetPr(p_OUT);
    aPr[0] = p.x;
    aPr[1] = p.y;

    return;
}
