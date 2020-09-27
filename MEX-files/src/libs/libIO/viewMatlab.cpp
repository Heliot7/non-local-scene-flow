#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "engine.h"

#include "libIO/viewMatlab.hpp"

using namespace std;

void showImage(unsigned char* img, int imgRows, int imgCols, int nChannels)
{
    Engine* ep;
    if (!(ep = engOpen("\0")))
    {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
    }
    else
    {
        const mwSize dims[3] = {imgRows, imgCols, nChannels};
        mxArray* M = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
        unsigned char* Mpr = (unsigned char*)mxGetData(M);
        memcpy(Mpr, img, sizeof(img)*imgRows*imgCols*nChannels);

        // Place variable into MATLAB workspace & plot the result
        engPutVariable(ep, "M", M);
        if(nChannels == 3)
        {
            engEvalString(ep, "figure; imshow(M);");
        }
        else if(nChannels == 1)
        {
            engEvalString(ep, "figure; imagesc(M);");
            engEvalString(ep, "colormap(gray);");
        }

        // Pause execution (does not work for QtCreator 1.3.1)
        unsigned char continueChar;
        cout << "Press any KEY to continue..." << endl;
        cin >> continueChar;

        // Free memory and close MATLAB
        mxDestroyArray(M);
        engEvalString(ep, "close;");
    }
}

void showMatrix(double* mat, int matRows, int matCols)
{
    Engine* ep;
    if (!(ep = engOpen("\0")))
    {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
    }
    else
    {
        mxArray* M = mxCreateDoubleMatrix(matRows, matCols, mxREAL);
        double* Mpr = mxGetPr(M);
        memcpy(Mpr, mat, sizeof(mat)*matRows*matCols);

        // Place variable into MATLAB workspace & plot the result
        engPutVariable(ep, "M", M);
        engEvalString(ep, "figure; imagesc(M);");
    }
}

