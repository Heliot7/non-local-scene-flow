#include <math.h>
#include <limits>
#include <cstring>
#include <algorithm>
#include <vector>

#include "mex.h" // for: mexPrintf(...)

#include "libIO/viewMatlab.hpp"

#include "libMath/matrix.hpp"
#include "libMath/metric.hpp"
#include "libImgProc/differentiation.hpp"
#include "libImgProc/filters.hpp"
#include "libImgProc/segmentation.hpp"

#include "nonLocalSceneFlow.hpp"

using namespace std;

void NonLocalSceneFlow::run()
{
    cout << "Hola" << endl;
}

// Wraps the run method within a MexFile
Result nonLocalSceneFlowMatlab(Input input)
{
    NonLocalSceneFlow* sceneFlow = new NonLocalSceneFlow(input);
    sceneFlow->run();
}

// Builds the input matrix of SOR method with equation values.
// INPUT PARAMETERS:
// - Set of warped images from a given view towards reference view at 't'
// - Set of warped images from a given view towards reference view at 't+1'
// - Set of partial deriviative_unknowns of all imageintensities at 't'
// - Set of partial deriviative_unknowns of all image intensities at 't+1'
// - Matrix of all fixed Psi robust coefficients for data term
// - Matrix of all fixed Psi diffusivity coefficients for smoothness term
// - Set of occlusion binary maps for all cameras w.r.t. reference image
void buildMatricesAb(double* mA, double* mB, int numCameras, int imgRows, int imgCols, int imgChannels, double* Z, double* u, double* v, double* w,
                             double* t0WarpedImgs, double* t1WarpedImgs, double* t0I_unknown, double* t1I_unknown, double* psiData, double* psiSmooth,
                             bool* occlusionBinary, double LEVEL, int numLevels, double factor, Weights weights, double* segmentPenalties)
{
    // Scaled values of alpha.
    double alphaZ = weights.alphaZ * pow(factor, (double)numLevels - LEVEL);
    double alphaUV = weights.alphaUV * pow(factor, (double)numLevels - LEVEL);
    double alphaW = weights.alphaW * pow(factor, (double)numLevels - LEVEL);

    // Initialise all matrices 'A' and 'b' for each pixel.
    for(int i = 0; i < imgRows*imgCols*4*4; ++i)
        mA[i] = 0.0;
    for(int i = 0; i < imgRows*imgCols*4; ++i)
        mB[i] = 0.0;

    // Limits...
    // 1 pxl due to central 2nd order differences.
    // int limit0 = 0;
    // int limitR = imgRows-1;
    // int limitC = imgCols-1;

    int offset = numCameras*imgRows*imgCols*imgChannels;
    int offsetOcc = numCameras*imgRows*imgCols;

    for(int COL = 1; COL < imgCols-1; ++COL)
        for(int ROW = 1 ; ROW < imgRows-1; ++ROW)
        {
            // Skip borders.
            //if(ROW <= limit0 || ROW >= limitR || COL <= limit0 || COL >= limitC)
            //    continue;

            int idx = COL*imgRows + ROW;

            // -> Segmentation info:
            double pUp = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*0];
            double pDown = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*1];
            double pLeft = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*2];
            double pRight = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*3];

            // First row of matrix (equation Z).
            Unknowns terms;
            terms.initZero();
            for(int CAMi = 0; CAMi < numCameras; ++CAMi)
                for(int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                {
                    int idxCam = CAMi + numCameras*ROW + numCameras*imgRows*COL;
                    int idxChannel = CAMi + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    int idxRef = 0 + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;

                    if(CAMi > 0)
                    {
                        // - Brightness Constancy at frame t
                        if(occlusionBinary[idxCam + offsetOcc*0])
                        {
                            terms.konst = terms.konst + psiData[idxChannel + offset*0] *
                                    (t0WarpedImgs[idxChannel] - t0WarpedImgs[idxRef]) *
                                    (t0I_unknown[idxChannel + offset*0] - t0I_unknown[idxRef + offset*0]);
                            terms.Z = terms.Z + psiData[idxChannel + offset*0] *
                                    (t0I_unknown[idxChannel + offset*0] - t0I_unknown[idxRef + offset*0]) *
                                    (t0I_unknown[idxChannel + offset*0] - t0I_unknown[idxRef + offset*0]);
                        }
                        // - Brightness Constancy at frame t+1
                        if(occlusionBinary[idxCam + offsetOcc*1])
                        {
                            terms.konst = terms.konst + psiData[idxChannel + offset*1] *
                                    (t1WarpedImgs[idxChannel] - t1WarpedImgs[idxRef]) *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]);
                            terms.Z = terms.Z + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]) *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]);
                            terms.u = terms.u + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]) *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]);
                            terms.v = terms.v + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]) *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]);
                            terms.w = terms.w + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]) *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]);
                        }
                    }
                    // - Brightness constancy between frames.
                    if(occlusionBinary[idxCam + offsetOcc*2])
                    {
                        terms.konst = terms.konst + psiData[idxChannel + offset*2] *
                                (t1WarpedImgs[idxChannel] - t0WarpedImgs[idxChannel]) *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]);
                        terms.Z = terms.Z + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]) *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]);
                        terms.u = terms.u + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*1]) *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]);
                        terms.v = terms.v + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*2]) *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]);
                        terms.w = terms.w + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*3]) *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]);
                    }
                }

            // - Add Smoothness term (divergence).

            // -> Diffusivity info:
            DiffNeighbours diff = computeDiffusivityZ(imgRows, imgCols, psiSmooth, ROW, COL);

            double divGrad = pUp*diff.up*(Z[COL*imgRows + ROW-1] - Z[COL*imgRows + ROW]) + pDown*diff.down*(Z[COL*imgRows + ROW+1] - Z[COL*imgRows + ROW]) +
                    pLeft*diff.left*(Z[(COL*imgRows-imgRows) + ROW] - Z[COL*imgRows + ROW]) + pRight*diff.right*(Z[(COL*imgRows+imgRows) + ROW] - Z[COL*imgRows + ROW]);
            double divIncr = pUp*diff.up + pDown*diff.down + pLeft*diff.left + pRight*diff.right;
            terms.konst = terms.konst - (alphaZ*weights.muZ/1.0)*(divGrad);
            terms.Z = terms.Z + (alphaZ*weights.muZ/1.0)*divIncr;

            // - Add all values into the matrices.
            mA[idx + imgRows*imgCols*0 + imgRows*imgCols*4*0] = terms.Z;
            mA[idx + imgRows*imgCols*0 + imgRows*imgCols*4*1] = terms.u;
            mA[idx + imgRows*imgCols*0 + imgRows*imgCols*4*2] = terms.v;
            mA[idx + imgRows*imgCols*0 + imgRows*imgCols*4*3] = terms.w;
            mB[idx + imgRows*imgCols*0] = -terms.konst;

            // 2nd row of matrix (equation u).
            terms.initZero();
            for(int CAMi = 0; CAMi < numCameras; ++CAMi)
                for(int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                {
                    int idxCam = CAMi + numCameras*ROW + numCameras*imgRows*COL;
                    int idxChannel = CAMi + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    int idxRef = 0 + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;

                    if(CAMi > 0)
                    {
                        // - Brightness Constancy: At frame 't+1'.
                        if(occlusionBinary[idxCam + offsetOcc*1])
                        {
                            terms.konst = terms.konst + psiData[idxChannel + offset*1] *
                                    (t1WarpedImgs[idxChannel] - t1WarpedImgs[idxRef]) *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]);
                            terms.Z = terms.Z + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]) *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]);
                            terms.u = terms.u + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]) *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]);
                            terms.v = terms.v + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]) *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]);
                            terms.w = terms.w + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]) *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]);
                        }
                    }
                    // - Brightness constancy between frames.
                    if(occlusionBinary[idxCam + offsetOcc*2])
                    {
                        terms.konst = terms.konst + psiData[idxChannel + offset*2] *
                                (t1WarpedImgs[idxChannel] - t0WarpedImgs[idxChannel]) *
                                (t1I_unknown[idxChannel + offset*1]);
                        terms.Z = terms.Z + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]) *
                                (t1I_unknown[idxChannel + offset*1]);
                        terms.u = terms.u + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*1]) *
                                (t1I_unknown[idxChannel + offset*1]);
                        terms.v = terms.v + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*2]) *
                                (t1I_unknown[idxChannel + offset*1]);
                        terms.w = terms.w + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*3]) *
                                (t1I_unknown[idxChannel + offset*1]);
                    }
                }

            // - Add Smoothness term (divergence).

            // -> Diffusivity info:
            diff = computeDiffusivityUVW(imgRows, imgCols, psiSmooth, ROW, COL);

            divGrad = pUp*diff.up*(u[COL*imgRows + ROW-1] - u[COL*imgRows + ROW]) + pDown*diff.down*(u[COL*imgRows + ROW+1] - u[COL*imgRows + ROW]) +
                    pLeft*diff.left*(u[(COL*imgRows-imgRows) + ROW] - u[COL*imgRows + ROW]) + pRight*diff.right*(u[(COL*imgRows+imgRows) + ROW] - u[COL*imgRows + ROW]);
            divIncr = pUp*diff.up + pDown*diff.down + pLeft*diff.left + pRight*diff.right;
            terms.konst = terms.konst - (alphaUV*weights.muUVW/1.0)*(divGrad);
            terms.u = terms.u + (alphaUV*weights.muUVW/1.0)*divIncr;

            // - Add all values into the matrices.
            mA[idx + imgRows*imgCols/**1 + imgRows*imgCols*4*0*/] = terms.Z;
            mA[idx + imgRows*imgCols/**1*/ + imgRows*imgCols*4*1] = terms.u;
            mA[idx + imgRows*imgCols/**1*/ + imgRows*imgCols*4*2] = terms.v;
            mA[idx + imgRows*imgCols/**1*/ + imgRows*imgCols*4*3] = terms.w;
            mB[idx + imgRows*imgCols/**1*/] = -terms.konst;

            // 3rd row of matrix (equation v).
            terms.initZero();
            for(int CAMi = 0; CAMi < numCameras; ++CAMi)
                for(int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                {
                    int idxCam = CAMi + numCameras*ROW + numCameras*imgRows*COL;
                    int idxChannel = CAMi + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    int idxRef = 0 + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;

                    if(CAMi > 0)
                    {
                        // - Brightness Constancy: At frame 't+1'.
                        if(occlusionBinary[idxCam + offsetOcc*1])
                        {
                            terms.konst = terms.konst + psiData[idxChannel + offset*1] *
                                    (t1WarpedImgs[idxChannel] - t1WarpedImgs[idxRef]) *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]);
                            terms.Z = terms.Z + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]) *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]);
                            terms.u = terms.u + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]) *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]);
                            terms.v = terms.v + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]) *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]);
                            terms.w = terms.w + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]) *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]);
                        }
                    }
                    // - Brightness constancy between frames.
                    if(occlusionBinary[idxCam + offsetOcc*2])
                    {
                        terms.konst = terms.konst + psiData[idxChannel + offset*2] *
                                (t1WarpedImgs[idxChannel] - t0WarpedImgs[idxChannel]) *
                                (t1I_unknown[idxChannel + offset*2]);
                        terms.Z = terms.Z + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]) *
                                (t1I_unknown[idxChannel + offset*2]);
                        terms.u = terms.u + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*1]) *
                                (t1I_unknown[idxChannel + offset*2]);
                        terms.v = terms.v + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*2]) *
                                (t1I_unknown[idxChannel + offset*2]);
                        terms.w = terms.w + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*3]) *
                                (t1I_unknown[idxChannel + offset*2]);
                    }
                }

            // -> Diffusivity info: (redundant)
            // diff = computeDiffusivityUVW(imgRows, imgCols, psiSmooth, ROW, COL);

            divGrad = pUp*diff.up*(v[COL*imgRows + ROW-1] - v[COL*imgRows + ROW]) + pDown*diff.down*(v[COL*imgRows + ROW+1] - v[COL*imgRows + ROW]) +
                    pLeft*diff.left*(v[(COL*imgRows-imgRows) + ROW] - v[COL*imgRows + ROW]) + pRight*diff.right*(v[(COL*imgRows+imgRows) + ROW] - v[COL*imgRows + ROW]);
            divIncr = pUp*diff.up + pDown*diff.down + pLeft*diff.left + pRight*diff.right;
            terms.konst = terms.konst - (alphaUV*weights.muUVW/1.0)*(divGrad);
            terms.v = terms.v + (alphaUV*weights.muUVW/1.0)*divIncr;

            // - Add all values into the matrices.
            mA[idx + imgRows*imgCols*2 + imgRows*imgCols*4*0] = terms.Z;
            mA[idx + imgRows*imgCols*2 + imgRows*imgCols*4*1] = terms.u;
            mA[idx + imgRows*imgCols*2 + imgRows*imgCols*4*2] = terms.v;
            mA[idx + imgRows*imgCols*2 + imgRows*imgCols*4*3] = terms.w;
            mB[idx + imgRows*imgCols*2] = -terms.konst;

            // 4th row of matrix (equation w).
            terms.initZero();
            for(int CAMi = 0; CAMi < numCameras; ++CAMi)
                for(int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                {
                    int idxCam = CAMi + numCameras*ROW + numCameras*imgRows*COL;
                    int idxChannel = CAMi + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    int idxRef = 0 + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;

                    if(CAMi > 0)
                    {
                        // - Brightness Constancy: At frame 't+1'.
                        if(occlusionBinary[idxCam + offsetOcc*1])
                        {
                            terms.konst = terms.konst + psiData[idxChannel + offset*1] *
                                    (t1WarpedImgs[idxChannel] - t1WarpedImgs[idxRef]) *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]);
                            terms.Z = terms.Z + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*0] - t1I_unknown[idxRef + offset*0]) *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]);
                            terms.u = terms.u + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*1] - t1I_unknown[idxRef + offset*1]) *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]);
                            terms.v = terms.v + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*2] - t1I_unknown[idxRef + offset*2]) *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]);
                            terms.w = terms.w + psiData[idxChannel + offset*1] *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]) *
                                    (t1I_unknown[idxChannel + offset*3] - t1I_unknown[idxRef + offset*3]);
                        }
                    }
                    // - Brightness constancy between frames.
                    if(occlusionBinary[idxCam + offsetOcc*2])
                    {
                        terms.konst = terms.konst + psiData[idxChannel + offset*2] *
                                (t1WarpedImgs[idxChannel] - t0WarpedImgs[idxChannel]) *
                                (t1I_unknown[idxChannel + offset*3]);
                        terms.Z = terms.Z + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*0] - t0I_unknown[idxChannel + offset*0]) *
                                (t1I_unknown[idxChannel + offset*3]);
                        terms.u = terms.u + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*1]) *
                                (t1I_unknown[idxChannel + offset*3]);
                        terms.v = terms.v + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*2]) *
                                (t1I_unknown[idxChannel + offset*3]);
                        terms.w = terms.w + psiData[idxChannel + offset*2] *
                                (t1I_unknown[idxChannel + offset*3]) *
                                (t1I_unknown[idxChannel + offset*3]);
                    }
                }

            // -> Diffusivity info: (redundant)
            // diff = computeDiffusivityUVW(imgRows, imgCols, psiSmooth, ROW, COL);

            divGrad = pUp*diff.up*(w[COL*imgRows + ROW-1] - w[COL*imgRows + ROW]) + pDown*diff.down*(w[COL*imgRows + ROW+1] - w[COL*imgRows + ROW]) +
                    pLeft*diff.left*(w[(COL*imgRows-imgRows) + ROW] - w[COL*imgRows + ROW]) + pRight*diff.right*(w[(COL*imgRows+imgRows) + ROW] - w[COL*imgRows + ROW]);
            divIncr = pUp*diff.up + pDown*diff.down + pLeft*diff.left + pRight*diff.right;
            terms.konst = terms.konst - (alphaW*weights.muUVW/1.0)*(divGrad);
            terms.w = terms.w + (alphaW*weights.muUVW/1.0)*divIncr;

            // - Add all values into the matrices.
            mA[idx + imgRows*imgCols*3 + imgRows*imgCols*4*0] = terms.Z;
            mA[idx + imgRows*imgCols*3 + imgRows*imgCols*4*1] = terms.u;
            mA[idx + imgRows*imgCols*3 + imgRows*imgCols*4*2] = terms.v;
            mA[idx + imgRows*imgCols*3 + imgRows*imgCols*4*3] = terms.w;
            mB[idx + imgRows*imgCols*3] = -terms.konst;
        }
}

// Private function: Calculates the diffusivity coefficient between a given pixel and all the
// 4-neighbours for unknown Z (3D reconstruction). This is used for the
// diffusivity term Psi'() in the smoothness term.
// INPUT PARAMETERS:
// - Matrix with all Diffusivity terms Psi'().
// - Pixel to work with.
DiffNeighbours computeDiffusivityZ(int imgRows, int imgCols, double* psiSmooth, int ROW, int COL)
{
    DiffNeighbours diff;
    diff.initZero();

    // TOP PIXEL:
    if(ROW != 1)
        diff.up = (psiSmooth[COL*imgRows + ROW] + psiSmooth[COL*imgRows + ROW-1]) / 2.0;
    else
        diff.up = (psiSmooth[COL*imgRows + ROW] + psiSmooth[COL*imgRows + ROW+1]) / 2.0;
    // DOWN PIXEL:
    if(ROW != imgRows-2)
        diff.down = (psiSmooth[COL*imgRows + ROW] + psiSmooth[COL*imgRows + ROW+1]) / 2.0;
    else
        diff.down = (psiSmooth[COL*imgRows + ROW] + psiSmooth[COL*imgRows + ROW-1]) / 2.0;
    // LEFT PIXEL:
    if(COL != 1)
        diff.left = (psiSmooth[COL*imgRows + ROW] + psiSmooth[(COL*imgRows-imgRows) + ROW]) / 2.0;
    else
        diff.left = (psiSmooth[COL*imgRows + ROW] + psiSmooth[(COL*imgRows+imgRows) + ROW]) / 2.0;
    // RIGHT PIXEL:
    if(COL != imgCols-2)
        diff.right = (psiSmooth[COL*imgRows + ROW] + psiSmooth[(COL*imgRows+imgRows) + ROW] ) / 2.0;
    else
        diff.right = (psiSmooth[COL*imgRows + ROW] + psiSmooth[(COL*imgRows-imgRows) + ROW] ) / 2.0;

    return diff;
}

// Private function: Calculates the diffusivity coefficient between a given pixel and all the
// 4-neighbours for disparity unknowns u, v and w. This is used for the
// diffusivity term Psi'() in the smoothness term.
// INPUT PARAMETERS:
// - Matrix with all Diffusivity terms Psi'().
// - Pixel to work with.
DiffNeighbours computeDiffusivityUVW(int imgRows, int imgCols, double* psiSmooth, int ROW, int COL)
{
    DiffNeighbours diff;
    diff.initZero();

    int offset = imgRows*imgCols;

    // TOP PIXEL:
    if(ROW != 1)
        diff.up = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[COL*imgRows + ROW-1 + offset]) / 2.0;
    else
        diff.up = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[COL*imgRows + ROW+1 + offset]) / 2.0;
    // DOWN PIXEL:
    if(ROW != imgRows-2)
        diff.down = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[COL*imgRows + ROW+1 + offset]) / 2.0;
    else
        diff.down = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[COL*imgRows + ROW-1 + offset]) / 2.0;
    // LEFT PIXEL:
    if(COL != 1)
        diff.left = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[(COL*imgRows-imgRows) + ROW + offset]) / 2.0;
    else
        diff.left = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[(COL*imgRows+imgRows) + ROW + offset]) / 2.0;
    // RIGHT PIXEL:
    if(COL != imgCols-2)
        diff.right = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[(COL*imgRows+imgRows) + ROW + offset] ) / 2.0;
    else
        diff.right = (psiSmooth[COL*imgRows + ROW + offset] + psiSmooth[(COL*imgRows-imgRows) + ROW + offset] ) / 2.0;
    
    return diff;
}

// Returns the Psi robust function derivative for the divergence part.
// INPUT PARAMETERS:
// - matrix with all partial derivates w.r.t. x of Z.
// - matrix with all partial derivates w.r.t. y of Z.
// - matrix with all partial derivates w.r.t. x of u.
// - matrix with all partial derivates w.r.t. y of u.
// - matrix with all partial derivates w.r.t. x of v.
// - matrix with all partial derivates w.r.t. y of v.
// - matrix with all partial derivates w.r.t. x of w.
// - matrix with all partial derivates w.r.t. y of w.
// - matrix with all partial derivates w.r.t. x of dZ.
// - matrix with all partial derivates w.r.t. y of dZ.
// - matrix with all partial derivates w.r.t. x of du.
// - matrix with all partial derivates w.r.t. y of du.
// - matrix with all partial derivates w.r.t. x of dv.
// - matrix with all partial derivates w.r.t. y of dv.
// - matrix with all partial derivates w.r.t. x of dw.
// - matrix with all partial derivates w.r.t. y of dw.
void buildDiffusivitySmooth(double* diff, int imgRows, int imgCols, double* Z_x, double* Z_y, double* u_x, double* u_y, double* v_x, double* v_y, double* w_x, double* w_y,
                                 double* dZ_x, double* dZ_y, double* du_x, double* du_y, double* dv_x, double* dv_y, double* dw_x, double* dw_y,
                                 double epsilon, double muZ, double muUVW)
{
    // Hardcoded to give previous results. We think it is not necessary to use muZ again inside PsiDeriv...
    muZ = 1.0;
    unsigned int idx, offset;

    for(unsigned int COL = 0; COL < imgCols; ++COL)
        for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
        {
            idx = COL*imgRows + ROW;
            offset = imgRows*imgCols;

            // With gradient of Z and dZ
            diff[idx /*+ offset*0*/] = psiDeriv( /*muZ**/(Z_x[idx] + dZ_x[idx])*(Z_x[idx] + dZ_x[idx]) + (Z_y[idx] + dZ_y[idx])*(Z_y[idx] + dZ_y[idx]), epsilon );

            // With gradients of flow (u,v,w) and their increments.
            diff[idx + offset/**1*/] = psiDeriv( /*muUVW**/
                                                 (u_x[idx] + du_x[idx])*(u_x[idx] + du_x[idx])  + (u_y[idx] + du_y[idx])*(u_y[idx] + du_y[idx])  +
                                                 /*muUVW**/
                                                 (v_x[idx] + dv_x[idx])*(v_x[idx] + dv_x[idx])  + (v_y[idx] + dv_y[idx])*(v_y[idx] + dv_y[idx])  +
                                                 /*muUVW**/
                                                 (w_x[idx] + dw_x[idx])*(w_x[idx] + dw_x[idx])  + (w_y[idx] + dw_y[idx])*(w_y[idx] + dw_y[idx]),
                                                 epsilon );
        }
}

// Calculates the robustness factor of data term with fixed increments.
// INPUT PARAMETERS:
// - Increment of Z
// - Increment of u
// - Increment of v
// - Increment of w
// - Set of warped images from a given view towards reference view at 't'
// - Set of warped images from a given view towards reference view at 't+1'
// - Set of partial deriviative_unknowns of all imageintensities at 't'
// - Set of partial deriviative_unknowns of all image intensities at 't+1'
void buildRobustnessData(double* psiData, int numCameras, int imgRows, int imgCols, int imgChannels, double* dZ, double* du, double* dv, double* dw,
                                          double* t0WarpedImgs, double* t1WarpedImgs, double* t0I_unknown, double* t1I_unknown, double epsilon)
{
    unsigned int idx, idxChannel, offset = numCameras*imgRows*imgCols*imgChannels;
    double termCAMt0, termCAMt1, termCAMref0, termCAMref1;
    double pxlDZ, pxlDU, pxlDV, pxlDW;

    for(unsigned int COL = 0; COL < imgCols; ++COL)
        for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
        {
            idx = COL*imgRows + ROW;
            pxlDZ = dZ[idx];
            pxlDU = du[idx];
            pxlDV = dv[idx];
            pxlDW = dw[idx];

            for(unsigned int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                for(unsigned int CAMi = 0; CAMi < numCameras; ++CAMi)
                {
                    idxChannel = CAMi + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;

                    termCAMt0 = t0WarpedImgs[idxChannel] + t0I_unknown[idxChannel /*+ offset*0*/] * pxlDZ;
                    termCAMt1 = t1WarpedImgs[idxChannel] + t1I_unknown[idxChannel /*+ offset*0*/] * pxlDZ +
                                      t1I_unknown[idxChannel + offset/**1*/] * pxlDU + t1I_unknown[idxChannel + offset*2] * pxlDV + t1I_unknown[idxChannel + offset*3] * pxlDW;

                    if(CAMi == 0)
                    {
                        // Brightness constancy at frame t
                        psiData[idxChannel /*+ offset*0*/] = 0.0;

                        // Brightness constancy at frame t+1
                        psiData[idxChannel + offset/**1*/] = 0.0;

                        // Keep info of first camera
                        termCAMref0 = termCAMt0;
                        termCAMref1 = termCAMt1;
                    }
                    else
                    {
                        // Brightness constancy at frame t
                        psiData[idxChannel /*+ offset*0*/] = psiDeriv( (termCAMt0 - termCAMref0)*(termCAMt0 - termCAMref0), epsilon );

                        // Brightness constancy at frame t+1
                        psiData[idxChannel + offset/**1*/] = psiDeriv( (termCAMt1 - termCAMref1)*(termCAMt1 - termCAMref1), epsilon );
                    }
    
                    // Brightness constancty between frames (temp)
                    psiData[idxChannel + offset*2] = psiDeriv( (termCAMt1 - termCAMt0)*(termCAMt1 - termCAMt0), epsilon );
                }
        }
}

// Private function: Calculates the derivative of robust function Psi.
// INPUT PARAMETERS:
// - Input of the function
double psiDeriv(double s, double epsilon)
{
    return 1 / (2 * sqrt(s + epsilon*epsilon));
}

// Calculates the partial derivatives w.r.t. unknowns for all pixels_i.
// INPUT PARAMETERS:
// - Projection matrices for all cameras
// - Focal length for all cameras
// - Principal point for all cameras
// - Depth image with Z unkonwn for all pixels w.r.t. reference view
// - Partial derivative w.r.t. x of Z
// - Partial derivative w.r.t. y of Z
// - Disparity in u-direction matrix w.r.t. reference view
// - Partial derivative w.r.t. x of u
// - Partial derivative w.r.t. y of u
// - Disparity in v-direction matrix w.r.t. reference view
// - Partial derivative w.r.t. x of v
// - Partial derivative w.r.t. y of v
// - Disparity in w-direction matrix w.r.t. reference view
// - Partial derivative w.r.t. x of w
// - Partial derivative w.r.t. y of w
// - Partial derivative w.r.t. x of warped images at frame 't'
// - Partial derivative w.r.t. y of warped images at frame 't'
// - frame whether at step 0 (t) or 1 (t+1)
void computePartialUnknowns(double* partialUnknowns, int numCameras, int imgRows, int imgCols, int imgChannels, double** projMatrices, double* focals, double* principals, double* camPos,
                                  double* Z, double* Z_x, double* Z_y, double* u, double* u_x, double* u_y, double* v, double* v_x, double* v_y, double* w, double* w_x, double* w_y,
                                  double* warpedImgs_x, double* warpedImgs_y, int numFrame)
{
    double focalX = focals[0];
    double focalY = focals[1];
    P2D focal; focal.x = focalX; focal.y = focalY;
    double principalX = principals[0];
    double principalY = principals[1];
    P2D cam; cam.x = camPos[0]; cam.y = camPos[1];
    unsigned int idx, idxChannel, offset = numCameras*imgRows*imgCols*imgChannels;
    double* M;

    for (unsigned int COL = 0; COL < imgCols; ++COL)
        for (unsigned int ROW = 0; ROW < imgRows; ++ROW)
        {
            idx = COL*imgRows + ROW;

            for (unsigned int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                for(unsigned int CAMi = 0; CAMi < numCameras; ++CAMi)
                {
                    M = projMatrices[CAMi];

                    idxChannel = CAMi + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    partialUnknowns[idxChannel /*+ offset*0*/] = 0.0; // dZ
                    partialUnknowns[idxChannel + offset/**1*/] = 0.0; // du
                    partialUnknowns[idxChannel + offset*2] = 0.0; // dv
                    partialUnknowns[idxChannel + offset*3] = 0.0; // dw

                    // +1 to get same results as Matlab enumeration from 1 till imgRows (not C++ like from 0 to imgRows-1)
                    P2D XYnoZ = getXYnoZ(COL+1, ROW+1, focalX, focalY, principalX, principalY);

                    P2D warpedImgs, t0I_iGradient_xy, t1I_iGradient_xy;
                    warpedImgs.x = warpedImgs_x[idxChannel];
                    warpedImgs.y = warpedImgs_y[idxChannel];

                    // Partial derivatives of pixels w.r.t. all 4 unknowns.
                    P2D pxl_iPartial_Z, pxl_iPartial_u, pxl_iPartial_v, pxl_iPartial_w;

                    if(numFrame == 0)
                    {
                        // Jacobian Matrix at frame 't' for current camera.
                        double J0[4];
                        buildJacobianT0(J0, M, focal, cam, XYnoZ.x, XYnoZ.y, Z[idx], Z_x[idx], Z_y[idx]);

                        // Gradient of a pixel w.r.t. reference (x,y) at 't'
                        t0I_iGradient_xy = p_iWarpedGradient(warpedImgs, J0);

                        // Partial derivative w.r.t. unknown of I_i (any camera).
                        pxl_iPartial_Z = t0PxlPartial_Z(M, cam, XYnoZ.x, XYnoZ.y, Z[idx]);
                        partialUnknowns[idxChannel + offset*0] = t0I_iGradient_xy.x * pxl_iPartial_Z.x + t0I_iGradient_xy.y * pxl_iPartial_Z.y;
                    }

                    else if(numFrame == 1)
                    {
                        // Jacobian Matrix at frame 't+1' for current camera.
                        double J1[4];
                        buildJacobianT1(J1, M, focal, cam, XYnoZ.x, XYnoZ.y, Z[idx], u[idx], v[idx], w[idx], Z_x[idx], Z_y[idx], u_x[idx], u_y[idx], v_x[idx], v_y[idx], w_x[idx], w_y[idx]);

                        // Gradient of a pixel w.r.t. reference (x,y) at 't+1'
                        t1I_iGradient_xy = p_iWarpedGradient(warpedImgs, J1);

                        // Partial derivative w.r.t. unknown of I_i (any camera).
                        pxl_iPartial_Z = t1PxlPartial_Z(M, cam, XYnoZ.x, XYnoZ.y, Z[idx], u[idx], v[idx], w[idx]);
                        partialUnknowns[idxChannel /*+ offset*0*/] = t1I_iGradient_xy.x * pxl_iPartial_Z.x + t1I_iGradient_xy.y * pxl_iPartial_Z.y;
                        pxl_iPartial_u = t1PxlPartial_u(M, cam, XYnoZ.x, XYnoZ.y, Z[idx], u[idx], v[idx], w[idx]);
                        partialUnknowns[idxChannel + offset/**1*/] = t1I_iGradient_xy.x * pxl_iPartial_u.x + t1I_iGradient_xy.y * pxl_iPartial_u.y;
                        pxl_iPartial_v = t1PxlPartial_v(M, cam, XYnoZ.x, XYnoZ.y, Z[idx], u[idx], v[idx], w[idx]);
                        partialUnknowns[idxChannel + offset*2] = t1I_iGradient_xy.x * pxl_iPartial_v.x + t1I_iGradient_xy.y * pxl_iPartial_v.y;
                        pxl_iPartial_w = t1PxlPartial_w(M, cam, XYnoZ.x, XYnoZ.y, Z[idx], u[idx], v[idx], w[idx]);
                        partialUnknowns[idxChannel + offset*3] = t1I_iGradient_xy.x * pxl_iPartial_w.x + t1I_iGradient_xy.y * pxl_iPartial_w.y;
                    }

                }
        }
}

// Private function: Builds a Jacobian matrix from a camera view 'i' to reference view at 't'.
// INPUT PARAMETERS:
// - projection matrix of camera 'i'
// - focal length of reference camera
// - value X without Z dependence
// - value Y without Z dependence
// - depth image with all Z values
// - partial derivative w.r.t. x of depth image with Z values
// - partial derivative w.r.t. y of depth image with Z values
void buildJacobianT0(double* J, double* M, P2D f, P2D cam, double X, double Y, double Z, double Z_x, double Z_y)
{
    // M (projMatrix) form:
    //   [0   3   6   9]
    //   [1   4   7  10]
    //   [2   5   8  11]

    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    double den = pow(MZ * Z + M[11], 2);

    // Partial derivative w.r.t x of x_i:
    double xx = ((M[0]/f.x)*Z + MX*Z_x)*(MZ*Z + M[11]) - (MX*Z + M[9] + M[0]*cam.x)*((M[2]/f.x)*Z + MZ*Z_x);
    double partial_xXi = xx / den;
    // Partial derivative w.r.t x of y_i:
    double xy = ((M[1]/f.x)*Z + MY*Z_x)*(MZ*Z + M[11]) - (MY*Z + M[10] + M[4]*cam.y)*((M[2]/f.x)*Z + MZ*Z_x);
    double partial_xYi = xy / den;
    // Partial derivative w.r.t y of x_i:
    double yx = ((M[3]/f.y)*Z + MX*Z_y)*(MZ*Z + M[11]) - (MX*Z + M[9] + M[0]*cam.x)*((M[2]/f.y)*Z + MZ*Z_y);
    double  partial_yXi = yx / den;
    // Partial derivative w.r.t y of y_i:
    double yy = ((M[4]/f.y)*Z + MY*Z_y)*(MZ*Z + M[11]) - (MY*Z + M[10] + M[4]*cam.y)*((M[2]/f.y)*Z + MZ*Z_y);
    double partial_yYi = yy / den;

    // Build 2x2 Jacobian matrix from partial derivatives:
    //   [0 partial_xXi 2 partial_yXi]
    //   [1 partial_xYi 3 partial_yYi]
    J[0] = partial_xXi;
    J[1] = partial_yXi;
    J[2] = partial_xYi;
    J[3] = partial_yYi;
}

// Private function: Builds a Jacobian matrix from a camera view 'i' to reference view at 't+1'.
// INPUT PARAMETERS:
// - current pixel
// - projection matrix of camera 'i'
// - focal length of reference camera
// - value X without Z dependence
// - value Y without Z dependence
// - depth image with all Z values
// - matrix with current u values
// - matrix with current v values
// - matrix with current w values
void buildJacobianT1(double* J, double* M, P2D f, P2D cam, double X, double Y, double Z, double u, double v, double w, double Z_x, double Z_y, double u_x, double u_y, double v_x, double v_y, double w_x, double w_y)
{
    // M (projMatrix) form:
    //   [0   3   6   9]
    //   [1   4   7  10]
    //   [2   5   8  11]

    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    double Mu = M[0] * u + M[3] * v + M[6] * w;
    double Mv = M[1] * u + M[4] * v + M[7] * w;
    double Mw = M[2] * u + M[5] * v + M[8] * w;
    double den = pow(MZ * Z + M[11], 2);

    // Partial derivative w.r.t x of x_i:
    double xx1 = ((M[0]/f.x)*Z + MX*Z_x + M[0]*u_x + M[3]*v_x + M[6]*w_x)*(MZ*Z + M[11] + Mw);
    double xx2 = (MX*Z + M[9] + M[0]*cam.x + Mu)*((M[2]/f.x)*Z + MZ*Z_x + M[2]*u_x + M[5]*v_x + M[8]*w_x);
    double partial_xXi = (xx1-xx2) / den;
    // Partial derivative w.r.t x of y_i:
    double xy1 = ((M[1]/f.x)*Z + MY*Z_x + M[1]*u_x + M[4]*v_x + M[7]*w_x)*(MZ*Z + M[11] + Mw);
    double xy2 = (MY*Z + M[10] + M[4]*cam.y + Mv)*((M[2]/f.x)*Z + MZ*Z_x + M[2]*u_x + M[5]*v_x + M[8]*w_x);
    double partial_xYi = (xy1-xy2) / den;
    // Partial derivative w.r.t y of x_i:
    double yx1 = ((M[3]/f.y)*Z + MX*Z_y + M[0]*u_y + M[3]*v_y + M[6]*w_y)*(MZ*Z + M[11] + Mw);
    double yx2 = (MX*Z + M[9] + M[0]*cam.x + Mu)*((M[2]/f.y)*Z + MZ*Z_y + M[2]*u_y + M[5]*v_y + M[8]*w_y);
    double partial_yXi = (yx1-yx2) / den;
    // Partial derivative w.r.t y of y_i:
    double yy1 = ((M[4]/f.y)*Z + MY*Z_y + M[1]*u_y + M[4]*v_y + M[7]*w_y)*(MZ*Z + M[11] + Mw);
    double yy2 = (MY*Z + M[10] + M[4]*cam.y + Mv)*((M[2]/f.y)*Z + MZ*Z_y + M[2]*u_y + M[5]*v_y + M[8]*w_y);
    double partial_yYi = (yy1-yy2) / den;

    // Build 2x2 Jacobian matrix from partial derivatives:
    //   [0 partial_xXi 2 partial_yXi]
    //   [1 partial_xYi 3 partial_yYi]
    J[0] = partial_xXi;
    J[1] = partial_yXi;
    J[2] = partial_xYi;
    J[3] = partial_yYi;
}

// Private function: Calculates the gradient at the given view but w.r.t. reference.
// INPUT PARAMETERS:
// - Gradient of warped image from given camera to reference.
// - Jacobian matrix to perform coordinate transformation.
P2D p_iWarpedGradient(P2D gradPxlRefImg, double* J)
{
    P2D gradient;
    if(fabs(J[0]*J[3] - J[2]*J[1]) < 1e-5)
    {
        gradient.x = 0.0;
        gradient.y = 0.0;
    }
    else
    {
        double invJ[4];
        inverse2x2(invJ, J);
        gradient.x = gradPxlRefImg.x * invJ[0] + gradPxlRefImg.y * invJ[2];
        gradient.y = gradPxlRefImg.x * invJ[1] + gradPxlRefImg.y * invJ[3];
    }
    return gradient;
}

// Private function: Returns partial derivative w.r.t. unknowns of a pixel_i.
// NOTE: Coming private functions differ in matters of unknown and frame.
// INPUT PARAMETERS:
// - Projection matrix of given view
// - 3D value at X axis
// - 3D value at Y axis
// - 3D value at Z axis
P2D t0PxlPartial_Z(double* M, P2D cam, double X, double Y, double Z)
{
    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    // Derivative calculation.
    P2D partial;
    double den = pow( MZ * Z + M[11], 2);
    double x_num = MX * (MZ * Z + M[11]) - (MX * Z + M[9] + M[0]*cam.x) * MZ;
    double y_num = MY * (MZ * Z + M[11]) - (MY * Z + M[10] + M[4]*cam.y) * MZ;
    partial.x = x_num / den;
    partial.y = y_num / den;

    return partial;
}
P2D t1PxlPartial_Z(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w)
{
    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    double Mu = M[0] * u + M[3] * v + M[6] * w;
    double Mv = M[1] * u + M[4] * v + M[7] * w;
    double Mw = M[2] * u + M[5] * v + M[8] * w;
    // Derivative calculation.
    P2D partial;
    double den = pow(MZ * Z + M[11] + Mw, 2);
    double x_num = MX * (MZ * Z + M[11] + Mw) - (MX * Z + M[9] + M[0]*cam.x + Mu) * MZ;
    double y_num = MY * (MZ * Z + M[11] + Mw) - (MY * Z + M[10] + M[4]*cam.y + Mv) * MZ;
    partial.x = x_num / den;
    partial.y = y_num / den;

    return partial;
}
P2D t1PxlPartial_u(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w)
{
    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    double Mu = M[0] * u + M[3] * v + M[6] * w;
    double Mv = M[1] * u + M[4] * v + M[7] * w;
    double Mw = M[2] * u + M[5] * v + M[8] * w;
    // Derivative calculation.
    P2D partial;
    double den = pow(MZ * Z + M[11] + Mw, 2);
    double x_num = M[0] * (MZ * Z + M[11] + Mw) - (MX * Z + M[9] + M[0]*cam.x + Mu) * M[2];
    double y_num = M[1] * (MZ * Z + M[11] + Mw) - (MY * Z + M[10] + M[4]*cam.y + Mv) * M[2];
    partial.x = x_num / den;
    partial.y = y_num / den;

    return partial;
}
P2D t1PxlPartial_v(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w)
{
    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    double Mu = M[0] * u + M[3] * v + M[6] * w;
    double Mv = M[1] * u + M[4] * v + M[7] * w;
    double Mw = M[2] * u + M[5] * v + M[8] * w;
    // Derivative calculation.
    P2D partial;
    double den = pow(MZ * Z + M[11] + Mw, 2);
    double x_num = M[3] * (MZ * Z + M[11] + Mw) - (MX * Z + M[9] + M[0]*cam.x + Mu) * M[5];
    double y_num = M[4] * (MZ * Z + M[11] + Mw) - (MY * Z + M[10] + M[4]*cam.y + Mv) * M[5];
    partial.x = x_num / den;
    partial.y = y_num / den;

    return partial;
}
P2D t1PxlPartial_w(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w)
{
    // Define some values.
    double MX = M[0] * X + M[3] * Y + M[6];
    double MY = M[1] * X + M[4] * Y + M[7];
    double MZ = M[2] * X + M[5] * Y + M[8];
    double Mu = M[0] * u + M[3] * v + M[6] * w;
    double Mv = M[1] * u + M[4] * v + M[7] * w;
    double Mw = M[2] * u + M[5] * v + M[8] * w;
    // Derivative calculation.
    P2D partial;
    double den = pow(MZ * Z + M[11] + Mw, 2);
    double x_num = M[6] * (MZ * Z + M[11] + Mw) - (MX * Z + M[9] + M[0]*cam.x + Mu) * M[8];
    double y_num = M[7] * (MZ * Z + M[11] + Mw) - (MY * Z + M[10] + M[4]*cam.y + Mv) * M[8];
    partial.x = x_num / den;
    partial.y = y_num / den;

    return partial;
}

// Private function: Calculates X' and Y', expressions without Z term.
// INPUT PARAMETERS:
// - Given pixel of current image
// - Scaled focal length
// - Principal point
P2D getXYnoZ(int posX, int posY, double focalX, double focalY, double principalX, double principalY)
{
    P2D XY;
    XY.x = posX/focalX - principalX/focalX;
    XY.y = posY/focalY - principalY/focalY;

    return XY;
}

// Private function: inverse of a matrix
void inverse2x2(double* invM, double* M)
{
    // Form of the 2x2 matrix
    // [ 0  2 ]
    // [ 1  3 ]
    double det = 1.0 / (M[0]*M[3] - M[2]*M[1]);
    invM[0] = M[3] * det;
    invM[1] = -M[1] * det;
    invM[2] = -M[2] * det;
    invM[3] = M[0] * det;
}

// Returns X and Y values of a 3D Point knowing camera parameters and Z.
// INPUT PARAMETERS:
// - Given pixel of current image
// - Current value of X.
// - Current value of Y.
// - Current value of Z.
// - Scaled focal length
// - Principal point
// - Camera position in world coordinates
void getXYfromZ(int imgRows, int imgCols, double* X, double* Y, double* Z, double* focals, double* principals, double* COP)
{
    int idx;

    for(unsigned int COL = 0; COL < imgCols; ++COL)
        for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
        {
            idx = COL*imgRows + ROW;

            X[idx] = Z[idx]*( (COL+1)/focals[0] - principals[0]/focals[0]) + COP[0];
            Y[idx] = Z[idx]*( (ROW+1)/focals[1] - principals[1]/focals[1]) + COP[1];
        }
}

// Calculates the projection of reference view into a given new view.
// INPUT PARAMETERS:
// - New view image at time 't'
// - Projection matrix of new view
// - 3D point w.r.t. each pixel
void buildProjectionMaps(double* projMap, int imgRows, int imgCols, int numCameras, double** projMatrices, double* X, double* Y, double* Z, double* u, double* v, double* w)
{
    double* projMatrix;
    int idx;
    P2D pxl;

    for(unsigned int COL = 0; COL < imgCols; ++COL)
        for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
        {
            idx = COL*imgRows + ROW;

            for(unsigned int CAM = 0; CAM < numCameras; ++CAM)
            {
                projMatrix = projMatrices[CAM];

                // Get pixel correspondence in the new view image (time t0 or t1)
                if(u == NULL)
                    pxl = projection3Dto2D(projMatrix, X[idx], Y[idx], Z[idx]);
                else
                    pxl = projection3Dto2D(projMatrix, X[idx]+u[idx], Y[idx]+v[idx], Z[idx]+w[idx]);

                // Store the pixel into the projection image correspondence.
                // NOTE: Due to row/col x/y mirroring, we swap results.
                projMap[CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*0] = pxl.y; // row
                projMap[CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*1] = pxl.x; // col
            }
        }
}

// Private function:
P2D projection3Dto2D(double* M, double X, double Y, double Z)
{
    P2D p;

    // p = (projMatrix(1:2,:)*[P;1]) / (projMatrix(3,:)*[P;1]);
    double den = (M[2]*X + M[5]*Y + M[8]*Z + M[11]);
    p.x = (M[0]*X + M[3]*Y + M[6]*Z + M[9])/den;
    p.y = (M[1]*X + M[4]*Y + M[7]*Z + M[10])/den;

    // Round at 4th decimal for precision after projection
    p.x = round(p.x*10000)/10000;
    p.y = round(p.y*10000)/10000;

    return p;
}

// Returns matrix with the occlusions of ref. view pixels given a new view.
// Possible occlusion map values at any index:
// 1 -> true, then adds new info to the energy functional.
// 0 -> false, therefore it does not add extra info to the energy functional.
// INPUT PARAMETERS:
// - Projection map of reference image pixels to corresponding camera view.
// - Matrix of reference image containing each pixel its point in 3D space.
// - Center of projection of corresponding camera view.
void buildOcclusionMaps(bool* occlusionMap, int numCameras, int imgRows, int imgCols, double* projImg, double* X, double* Y, double* Z, double* u, double* v, double* w, double** COP)
{
    // Initialise the future occlusion map matrix with all pixels occluded.
    for(unsigned int PXLi = 0; PXLi < numCameras*imgRows*imgCols; ++PXLi)
        occlusionMap[PXLi] = false;

    // Matrix with non-occluded ref. pixels in corresponding camera view.
    vector< vector< P2D > > visiblePxls;
    visiblePxls.resize(imgRows);
    for (int ROW = 0; ROW < imgRows; ++ROW)
        visiblePxls[ROW].resize(imgCols);

    for(int CAM = 0; CAM < numCameras; ++CAM)
    {
        // Get COP of corresponding camera
        double* camCOP = COP[CAM];

        // Initialise to [0 0].
        for(unsigned int COL = 0; COL < imgCols; ++COL)
            for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
            {
                P2D pxl;
                pxl.initZero();
                visiblePxls[ROW][COL] = pxl;
            }

        for(unsigned int COL = 0; COL < imgCols; ++COL)
            for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
            {
                int idx = COL*imgRows + ROW;

                // Get mapped pixel in corresponding camera view.
                P2D pxlNewCam;
                pxlNewCam.row = projImg[CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*0];
                pxlNewCam.col = projImg[CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*1];

                // To be visible, the point must lie within the image view.
                if(!isPxlOutOfBounds(imgRows, imgCols, pxlNewCam))
                {
                    P2D currVisible = visiblePxls[(int)round(pxlNewCam.row)-1][(int)round(pxlNewCam.col)-1];

                    // No previous pixels correspondences
                    if(currVisible.row == 0.0 && currVisible.col == 0.0)
                    {
                        P2D pxl; pxl.x = ROW+1; pxl.y = COL+1;
                        visiblePxls[(int)round(pxlNewCam.row)-1][(int)round(pxlNewCam.col)-1] = pxl;
                        occlusionMap[CAM + numCameras*ROW + numCameras*imgRows*COL] = true;
                    }
                    // There are other candidates, need to decide which stays.
                    else
                    {
                        // Get 3D point of current visible pixel.
                        int idxCurrent = imgRows*(int)currVisible.col + (int)currVisible.row;
                        double distCurr3D = sqrt((X[idxCurrent] - camCOP[0])*(X[idxCurrent] - camCOP[0]) + (Y[idxCurrent] - camCOP[1])*(Y[idxCurrent] - camCOP[1])
                                + (Z[idxCurrent] - camCOP[2])*(Z[idxCurrent] - camCOP[2]));
                        // Get 3D point of new candidate pixel.
                        double distNew3D = sqrt((X[idx] - camCOP[0])*(X[idx] - camCOP[0]) + (Y[idx] - camCOP[1])*(Y[idx] - camCOP[1]) + (Z[idx] - camCOP[2])*(Z[idx] - camCOP[2]));

                        // Calculate distance for each point and the COP.
                        // Smallest value is visible in the image (~occluded).
                        if(distNew3D < distCurr3D)
                        {
                            // New candidate is currently not occluded.
                            occlusionMap[CAM + numCameras*ROW + numCameras*imgRows*COL] = true;
                            // While old one is from now on occluded.
                            occlusionMap[CAM + numCameras*(int)currVisible.row + numCameras*imgRows*(int)currVisible.col] = false;
                            // We also need to update the visible pxl matrix.
                            P2D pxl; pxl.x = ROW+1; pxl.y = COL+1;
                            visiblePxls[(int)round(pxlNewCam.row)-1][(int)round(pxlNewCam.col)-1] = pxl;
                        }
                    }
                }
            }
    }
}

// Returns whether a pixel is out of an image region.
// INPUT PARAMETERS:
// - Image number of rows and columns
// - Pixel coordinate candidate.
bool isPxlOutOfBounds(int imgRows, int imgCols, P2D pxl)
{
    return (pxl.row < 1 || pxl.col < 1 || pxl.row > imgRows || pxl.col > imgCols);
}

// Returns a set of occlusionBinary maps for each views and equation part.
// INPUT PARAMETERS:
// - Occlusion map for images at frame 't'.
// - Occlusion map for images at frame 't+1'.
void buildOcclusionBinary(bool* occBin, int numCameras, int imgRows, int imgCols, double* t0ProjImg, double* t1ProjImg, bool* t0occlusionMaps, bool* t1occlusionMaps)
{
    unsigned idx, idxRef, offset = numCameras*imgRows*imgCols;
    P2D t0Pxl, t1Pxl;

    for(unsigned int COL = 0; COL < imgCols; ++COL)
        for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
            for(unsigned int CAM = 0; CAM < numCameras; ++CAM)
            {
                idx = CAM + numCameras*ROW + numCameras*imgRows*COL;
                idxRef = /*0 +*/ numCameras*ROW + numCameras*imgRows*COL; // Camera 0

                t0Pxl.row = t0ProjImg[idx /*+ offset*0*/];
                t0Pxl.col = t0ProjImg[idx + offset/**1*/];
                t1Pxl.row = t1ProjImg[idx /*+ offset*0*/];
                t1Pxl.col = t1ProjImg[idx + offset/**1*/];

                // "t0, t1 and temp for each CAM = TRUE" will not be treated in the future.
                occBin[idx /*+ offset*0*/] = t0occlusionMaps[idx] && t0occlusionMaps[idxRef] && !isPxlOutOfBounds(imgRows, imgCols, t0Pxl);
                occBin[idx + offset/**1*/] = t1occlusionMaps[idx] && t1occlusionMaps[idxRef] && !isPxlOutOfBounds(imgRows, imgCols, t1Pxl);
                occBin[idx + offset*2] = t0occlusionMaps[idx] && t1occlusionMaps[idx] && !isPxlOutOfBounds(imgRows, imgCols, t0Pxl) && !isPxlOutOfBounds(imgRows, imgCols, t1Pxl);
            }
}

// Returns a warped image in reference view from a give 'i' view.
// INPUT PARAMETERS:
// - Projection map in reference image pixels to corresponding camera view.
// - Image of corresponding camera view.
void buildWarpedImage(double* warpedImgs, int numCameras, int imgRows, int imgCols, int imgChannels, unsigned char** images, double* projImgs, bool* occlusionMaps)
{
    // Init to 0
    for(int PXLi = 0; PXLi < numCameras*imgRows*imgCols*imgChannels; ++PXLi)
        warpedImgs[PXLi] = 0.0;

    unsigned int idx, idxOcc, idxProjRow, idxProjCol, pROW, pCOL;
    double pxlProjRow, pxlProjCol;

    // Add intensity value from image to warped image at [ROW,COL].
    for(unsigned int CAM = 0; CAM < numCameras; ++CAM)
    {
        // Step 1: Warping of pixel position with a valid assignment
        for(unsigned int COL = 0; COL < imgCols; ++COL)
            for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
            {
                idxOcc = CAM + numCameras*ROW + numCameras*imgRows*COL;
                idxProjRow = CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*0;
                idxProjCol = CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*1;
                for(unsigned int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                {
                    idx = CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    if(occlusionMaps[idxOcc] == true)
                    {
                        pxlProjRow = projImgs[idxProjRow];
                        pxlProjCol = projImgs[idxProjCol];
                        pROW = (unsigned int)floor(pxlProjRow);
                        pCOL = (unsigned int)floor(pxlProjCol);
                        P2D p11, p12, p21, p22;
                        p11.x = pROW+1; p11.y = pCOL;
                        double i11 = (double)images[CAM][pCOL*imgRows-imgRows + pROW + CHANNEL*(imgRows*imgCols)];
                        p12.x = pROW; p12.y = pCOL;
                        double i12 = (double)images[CAM][pCOL*imgRows-imgRows + pROW-1 + CHANNEL*(imgRows*imgCols)];
                        p21.x = pROW+1; p21.y = pCOL+1;
                        double i21 = (double)images[CAM][pCOL*imgRows + pROW + CHANNEL*(imgRows*imgCols)];
                        p22.x = pROW; p22.y = pCOL+1;
                        double i22 = (double)images[CAM][pCOL*imgRows + pROW-1 + CHANNEL*(imgRows*imgCols)];
                        warpedImgs[idx] = bilinearInterpolation(pxlProjRow, pxlProjCol, p11, i11, p12, i12, p21, i21, p22, i22);
                    }
                }
            }

        // Step 2: If no projection from ref-new view, bilinear interpolation to fill it.
        for(unsigned int COL = 0; COL < imgCols; ++COL)
            for(unsigned int ROW = 0; ROW < imgRows; ++ROW)
            {
                idxOcc = CAM + numCameras*ROW + numCameras*imgRows*COL;
                for(unsigned int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
                {
                    idx = CAM + numCameras*ROW + numCameras*imgRows*COL + numCameras*imgRows*imgCols*CHANNEL;
                    if(occlusionMaps[idxOcc] == false && ROW > 0 && ROW < imgRows-1 && COL > 0 && COL < imgCols-1)
                    {
                        P2D p11, p12, p21, p22;
                        p11.x = ROW+2; p11.y = COL;
                        double i11 = warpedImgs[CAM + numCameras*ROW+numCameras + numCameras*imgRows*COL-numCameras*imgRows + numCameras*imgRows*imgCols*CHANNEL];
                        p12.x = ROW; p12.y = COL;
                        double i12 = warpedImgs[CAM + numCameras*ROW-numCameras + numCameras*imgRows*COL-numCameras*imgRows + numCameras*imgRows*imgCols*CHANNEL];
                        p21.x = ROW+2; p21.y = COL+2;
                        double i21 = warpedImgs[CAM + numCameras*ROW+numCameras + numCameras*imgRows*COL+numCameras*imgRows + numCameras*imgRows*imgCols*CHANNEL];
                        p22.x = ROW; p22.y = COL+2;
                        double i22 = warpedImgs[CAM + numCameras*ROW-numCameras + numCameras*imgRows*COL+numCameras*imgRows + numCameras*imgRows*imgCols*CHANNEL];
                        warpedImgs[idx] = bilinearInterpolation(ROW+1, COL+1, p11, i11, p12, i12, p21, i21, p22, i22);
                    }
                }
            }
    }
}

// Private function: Returns an intensity value calculated from bilinear interpolation
// INPUT PARAMETERS:
// - Bottom-left intensity value
// - Top-left intensity value
// - Bottom-right intensity value
// - Top-right intensity value

//  p12 ========== p22
//   ||             ||
//   ||             ||
//   ||      p      ||
//   ||             ||
//   ||             ||
//   p11 ========== p21
double bilinearInterpolation(double ROW, double COL, P2D p11, double i11, P2D p12, double i12, P2D p21, double i21, P2D p22, double i22)
{
    double x = COL;
    double y = ROW;
    double x1 = p11.y;
    double x2 = p22.y;
    double y1 = p21.x;
    double y2 = p12.x;
    // Step 1: X AXIS
    double R1 = ( (x2 - x) / (x2 - x1) ) * i11 + ( (x - x1) / (x2 - x1) ) * i21;
    double R2 = ( (x2 - x) / (x2 - x1) ) * i12 + ( (x - x1) / (x2 - x1) ) * i22;
    // Step 2: Y AXIS
    return ( (y2 - y) / (y2 - y1) ) * R1 + ( (y - y1) / (y2 - y1) ) * R2;
}

P2D projection3Dto2D(double* projMatrix, double* P)
{
    P2D p;

    double denom = projMatrix[2]*P[0] + projMatrix[5]*P[1] + projMatrix[8]*P[2] + projMatrix[11]*P[3];
    p.x = (projMatrix[0]*P[0] + projMatrix[3]*P[1] + projMatrix[6]*P[2] + projMatrix[9]*P[3]) / denom;
    p.y = (projMatrix[1]*P[0] + projMatrix[4]*P[1] + projMatrix[7]*P[2] + projMatrix[10]*P[3]) / denom;

    // Round at 4th decimal for precision after projection
    p.x = round(p.x*10000.0)/10000.0;
    p.y = round(p.y*10000.0)/10000.0;

    return p;
}

/* Iterative Solvers */

// Iteratives over all pixels to estimate increments with an iterative solver method.
// INPUT PARAMETERS:
// - Matrix of all A matrices per pixel
// - Matrix of all b matrices per pixel
// - Matrix of all diffusivity coefficients per pixel
// - Current level in the multiresolution pyramid
void iterativeSolver(double* dZ, double* du, double* dv, double* dw, double* evolResidual, char* itSolver, char* itPrecon, int imgRows, int imgCols, double* a, double* b, double* diffSmooth,
                                   double wSOR, int maxIter, double errorSolver, int LEVEL, int numLevels, double factor, Weights weights, double* segmentPenalties)
{
    double incr[4];

    // Declare stencils of matrix A
    double** aStencil = new double*[imgRows*imgCols];

    // -> 2-norm of all elems in matrix b.
    double valueB = 0.0;

    // SOR method:
    double bResidual[4], div_d[4];

    // CG method:
    double stencil[20], **resCG, **resPCG, **dirCG, **Adir, **Ares;
    // Preconditioning:
    double** M; double** preAdir;
    if(strcmp(itSolver, "CG") == 0 || strcmp(itSolver, "CR") == 0)
    {
        resCG = new double*[imgRows*imgCols];
        resPCG = new double*[imgRows*imgCols];
        dirCG = new double*[imgRows*imgCols];
        Adir = new double*[imgRows*imgCols];
        Ares = new double*[imgRows*imgCols];

        // Initialise to 0 all at once
        for(int ROW = 0; ROW < imgRows; ++ROW)
            for(int COL = 0; COL < imgCols; ++COL)
            {
                int idx = COL*imgRows + ROW;
                resCG[idx] = new double[4];
                resPCG[idx] = new double[4];
                Ares[idx] = new double[4];
                dirCG[idx] = new double[4];
                Adir[idx] = new double[4];
            }

        if(strcmp(itPrecon, "NO") != 0)
        {
            M = new double*[imgRows*imgCols];
            if(strcmp(itSolver, "CR") == 0)
                preAdir = new double*[imgRows*imgCols];
            for(int ROW = 0; ROW < imgRows; ++ROW)
                for(int COL = 0; COL < imgCols; ++COL)
                {
                    int idx = imgRows*COL + ROW;
                    M[idx] = new double[4*20];
                }

            if(strcmp(itSolver, "CR") == 0)
            {
                for(int ROW = 0; ROW < imgRows; ++ROW)
                    for(int COL = 0; COL < imgCols; ++COL)
                    {
                        int idx = imgRows*COL + ROW;
                        preAdir[idx] = new double[4];
                    }
            }
        }
    }

    // Scaled values of alpha.
    double alphaZ = weights.alphaZ * pow(factor, (double)numLevels - LEVEL);
    double alphaUV = weights.alphaUV * pow(factor, (double)numLevels - LEVEL);
    double alphaW = weights.alphaW * pow(factor, (double)numLevels - LEVEL);

    // State border limits due to finite differences
    int limit = 1;
    int rowPxl = imgRows*imgCols;
    int colPxl = imgRows*imgCols*4;

    // Declare stencil A matrix.
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = COL*imgRows + ROW;

            if(ROW < limit || COL < limit || ROW >= imgRows-limit || COL >= imgCols-limit)
            {
                if(strcmp(itSolver, "CG") == 0 || strcmp(itSolver, "CR") == 0)
                {
                    dirCG[idx][0] = 0.0; dirCG[idx][1] = 0.0; dirCG[idx][2] = 0.0; dirCG[idx][3] = 0.0;
                    resCG[idx][0] = 0.0; resCG[idx][1] = 0.0; resCG[idx][2] = 0.0; resCG[idx][3] = 0.0;
                }
                continue;
            }

            int factorBorders = 1;

            // -> Diffusivity info:
            DiffNeighbours diffZ = computeDiffusivityZ(imgRows, imgCols, diffSmooth, ROW, COL);
            DiffNeighbours diffUVW = computeDiffusivityUVW(imgRows, imgCols, diffSmooth, ROW, COL);

            // -> Segmentation info:
            double pUp = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*0];
            double pDown = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*1];
            double pLeft = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*2];
            double pRight = segmentPenalties[imgRows*COL + ROW + imgRows*imgCols*3];

            // Stencil matrix A [4,4x5]:
            aStencil[idx] = new double[4*20];
            for(int coef = 0; coef < 4*20; ++coef)
                aStencil[idx][coef] = 0.0;

            // -> middle cols: center pixel
            aStencil[idx][8*4+0] = a[idx + rowPxl*0 + colPxl*0]; aStencil[idx][9*4] = a[idx + rowPxl*0 + colPxl*1]; aStencil[idx][10*4] = a[idx + rowPxl*0 + colPxl*2]; aStencil[idx][11*4] = a[idx + rowPxl*0 + colPxl*3];
            aStencil[idx][8*4+1] = a[idx + rowPxl*1 + colPxl*0]; aStencil[idx][9*4+1] = a[idx + rowPxl*1 + colPxl*1]; aStencil[idx][10*4+1] = a[idx + rowPxl*1 + colPxl*2]; aStencil[idx][11*4+1] = a[idx + rowPxl*1 + colPxl*3];
            aStencil[idx][8*4+2] = a[idx + rowPxl*2 + colPxl*0]; aStencil[idx][9*4+2] = a[idx + rowPxl*2 + colPxl*1]; aStencil[idx][10*4+2] = a[idx + rowPxl*2 + colPxl*2]; aStencil[idx][11*4+2] = a[idx + rowPxl*2 + colPxl*3];
            aStencil[idx][8*4+3] = a[idx + rowPxl*3 + colPxl*0]; aStencil[idx][9*4+3] = a[idx + rowPxl*3 + colPxl*1]; aStencil[idx][10*4+3] = a[idx + rowPxl*3 + colPxl*2]; aStencil[idx][11*4+3] = a[idx + rowPxl*3 + colPxl*3];

            // -> 1st cols: top pixel
            if(ROW > limit)
            {
                aStencil[idx][0*4+0] -= pUp*(alphaZ*weights.muZ/1.0)*diffZ.up;
                aStencil[idx][1*4+1] -= pUp*(alphaUV*weights.muUVW/1.0)*diffUVW.up;
                aStencil[idx][2*4+2] -= pUp*(alphaUV*weights.muUVW/1.0)*diffUVW.up;
                aStencil[idx][3*4+3] -= pUp*(alphaW*weights.muUVW/1.0)*diffUVW.up;
            }
            else
            {
                factorBorders *= 2;
                aStencil[idx][16*4+0] -= pDown*(alphaZ*weights.muZ/1.0)*diffZ.down;
                aStencil[idx][17*4+1] -= pDown*(alphaUV*weights.muUVW/1.0)*diffUVW.down;
                aStencil[idx][18*4+2] -= pDown*(alphaUV*weights.muUVW/1.0)*diffUVW.down;
                aStencil[idx][19*4+3] -= pDown*(alphaW*weights.muUVW/1.0)*diffUVW.down;
            }
            // -> 2nd cols: left pixel
            if(COL > limit)
            {
                aStencil[idx][4*4+0] -= pLeft*(alphaZ*weights.muZ/1.0)*diffZ.left;
                aStencil[idx][5*4+1] -= pLeft*(alphaUV*weights.muUVW/1.0)*diffUVW.left;
                aStencil[idx][6*4+2] -= pLeft*(alphaUV*weights.muUVW/1.0)*diffUVW.left;
                aStencil[idx][7*4+3] -= pLeft*(alphaW*weights.muUVW/1.0)*diffUVW.left;
            }
            else
            {
                factorBorders *= 2;
                aStencil[idx][12*4+0] -= pRight*(alphaZ*weights.muZ/1.0)*diffZ.right;
                aStencil[idx][13*4+1] -= pRight*(alphaUV*weights.muUVW/1.0)*diffUVW.right;
                aStencil[idx][14*4+2] -= pRight*(alphaUV*weights.muUVW/1.0)*diffUVW.right;
                aStencil[idx][15*4+3] -= pRight*(alphaW*weights.muUVW/1.0)*diffUVW.right;
            }
            // -> 4th cols: right pixel
            if(COL < imgCols - limit - 1)
            {
                aStencil[idx][12*4+0] -= pRight*(alphaZ*weights.muZ/1.0)*diffZ.right;
                aStencil[idx][13*4+1] -= pRight*(alphaUV*weights.muUVW/1.0)*diffUVW.right;
                aStencil[idx][14*4+2] -= pRight*(alphaUV*weights.muUVW/1.0)*diffUVW.right;
                aStencil[idx][15*4+3] -= pRight*(alphaW*weights.muUVW/1.0)*diffUVW.right;
            }
            else
            {
                factorBorders *= 2;
                aStencil[idx][4*4+0] -= pLeft*(alphaZ*weights.muZ/1.0)*diffZ.left;
                aStencil[idx][5*4+1] -= pLeft*(alphaUV*weights.muUVW/1.0)*diffUVW.left;
                aStencil[idx][6*4+2] -= pLeft*(alphaUV*weights.muUVW/1.0)*diffUVW.left;
                aStencil[idx][7*4+3] -= pLeft*(alphaW*weights.muUVW/1.0)*diffUVW.left;
            }
            // -> 5th cols: down pixel
            if(ROW < imgRows - limit - 1)
            {
                aStencil[idx][16*4+0] -= pDown*(alphaZ*weights.muZ/1.0)*diffZ.down;
                aStencil[idx][17*4+1] -= pDown*(alphaUV*weights.muUVW/1.0)*diffUVW.down;
                aStencil[idx][18*4+2] -= pDown*(alphaUV*weights.muUVW/1.0)*diffUVW.down;
                aStencil[idx][19*4+3] -= pDown*(alphaW*weights.muUVW/1.0)*diffUVW.down;
            }
            else
            {
                factorBorders *= 2;
                aStencil[idx][0*4+0] -= pUp*(alphaZ*weights.muZ/1.0)*diffZ.up;
                aStencil[idx][1*4+1] -= pUp*(alphaUV*weights.muUVW/1.0)*diffUVW.up;
                aStencil[idx][2*4+2] -= pUp*(alphaUV*weights.muUVW/1.0)*diffUVW.up;
                aStencil[idx][3*4+3] -= pUp*(alphaW*weights.muUVW/1.0)*diffUVW.up;
            }

            // Adapt aStencil and b matrix in the border factor to be symmetric.
            for(int j = 0; j < 4*20; ++j)
                aStencil[idx][j] /= factorBorders;
            for(int j = 0; j < 4; ++j)
            {
                b[idx + rowPxl*j] /= factorBorders;
                valueB += b[idx + rowPxl*j]*b[idx + rowPxl*j];
            }

            if(strcmp(itSolver, "CG") == 0 || strcmp(itSolver, "CR") == 0)
            {
                // As incr = 0 always at 1st iter, r = b - Ax = b (as x = 0^T);
                resCG[idx][0] = b[idx + rowPxl*0]; resCG[idx][1] = b[idx + rowPxl*1]; resCG[idx][2] = b[idx + rowPxl*2]; resCG[idx][3] = b[idx + rowPxl*3];
                resPCG[idx][0] = resCG[idx][0]; resPCG[idx][1] = resCG[idx][1]; resPCG[idx][2] = resCG[idx][2]; resPCG[idx][3] = resCG[idx][3];
                if(strcmp(itSolver, "CR") == 0)
                {
                    resCG[idx][0] = resPCG[idx][0]; resCG[idx][1] = resPCG[idx][1]; resCG[idx][2] = resPCG[idx][2]; resCG[idx][3] = resPCG[idx][3];
                    Ares[idx][0] = resPCG[idx][0]; Ares[idx][1] = resPCG[idx][1]; Ares[idx][2] = resPCG[idx][2]; Ares[idx][3] = resPCG[idx][3];
                }
                dirCG[idx][0] = resPCG[idx][0]; dirCG[idx][1] = resPCG[idx][1]; dirCG[idx][2] = resPCG[idx][2]; dirCG[idx][3] = resPCG[idx][3];
                Adir[idx][0] = resPCG[idx][0]; Adir[idx][1] = resPCG[idx][1]; Adir[idx][2] = resPCG[idx][2]; Adir[idx][3] = resPCG[idx][3];
            }
        }

    // Apply preconditioner to residual
    if(strcmp(itPrecon, "NO") != 0)
    {
        // Create preconditioning matrix M:
        computePreconditioner(M, imgRows, imgCols, limit, itPrecon, aStencil);
        applyPreconditioner(imgRows, imgCols, limit, resPCG, resCG, M);

        for(int ROW = limit; ROW < imgRows-limit; ++ROW)
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int j = 0; j < 4; ++j)
                {
                    int i = imgRows*COL + ROW;
                    dirCG[i][j] = resPCG[i][j];
                    Adir[i][j] = resPCG[i][j];
                }

        if(strcmp(itSolver, "CR") == 0)
            for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                for(int COL = limit; COL < imgCols-limit; ++COL)
                    for(int j = 0; j < 4; ++j)
                    {
                        int i = imgRows*COL + ROW;
                        resCG[i][j] = resPCG[i][j];
                        Ares[i][j] = resPCG[i][j];
                    }
    }

    valueB = sqrt(valueB);

    // SOR
    int scanMode = 0;

    // Initialise residual for each pixel
    double** residual = new double*[imgRows*imgCols];
    for(int ROW = limit; ROW < imgRows-limit; ++ROW)
        for(int COL = limit; COL < imgCols-limit; ++COL)
        {
            int idx = imgRows*COL + ROW;
            residual[idx] = new double[4];
            residual[idx][0] = numeric_limits<double>::infinity();
            residual[idx][1] = numeric_limits<double>::infinity();
            residual[idx][2] = numeric_limits<double>::infinity();
            residual[idx][3] = numeric_limits<double>::infinity();
        }

    // Iterative a max number of times (or till convergence)
    int it = 0;
    bool convergence = false;
    while(!convergence && it < maxIter)
    {
        // mexPrintf("it: %d \n", it);

        // ==========================
        // Conjugate Gradient Method:
        // ==========================
        if(strcmp(itSolver, "CG") == 0)
        {
            // Updates multiplication A*p per pixel 4 unknowns
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;

                    for(int i = 0; i < 4; ++i)
                    {
                        // UP
                        stencil[i] = dirCG[COL*imgRows + (ROW-1)][i];
                        // LEFT
                        stencil[4 + i] = dirCG[(COL-1)*imgRows + ROW][i];
                        // CENTER
                        stencil[8 + i] = dirCG[COL*imgRows + ROW][i];
                        // RIGHT
                        stencil[12 + i] = dirCG[(COL+1)*imgRows + ROW][i];
                        // BOTTOM
                        stencil[16 + i] = dirCG[COL*imgRows + (ROW+1)][i];
                    }
                    mulAp(Adir[idx], aStencil[idx], stencil);
                }

            // Calculate Alpha value
            double numAlpha = 0.0;
            double denAlpha = 0.0;
            for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                for(int COL = limit; COL < imgCols-limit; ++COL)
                {
                    int idx = COL*imgRows + ROW;
                    numAlpha += resCG[idx][0]*resPCG[idx][0] + resCG[idx][1]*resPCG[idx][1] + resCG[idx][2]*resPCG[idx][2] + resCG[idx][3]*resPCG[idx][3];
                    denAlpha += dirCG[idx][0]*Adir[idx][0] + dirCG[idx][1]*Adir[idx][1] + dirCG[idx][2]*Adir[idx][2] + dirCG[idx][3]*Adir[idx][3];
                }
            double denBeta = numAlpha;
            double alpha = numAlpha / denAlpha;

            // Calculation of updated increments and residual
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    incr[0] = dZ[idx];
                    incr[1] = du[idx];
                    incr[2] = dv[idx];
                    incr[3] = dw[idx];

                    CG_incr_res(incr, resCG[idx], dirCG[idx], Adir[idx], alpha);
                    for(int i = 0; i < 4; ++i)
                        residual[idx][i] = resCG[idx][i];

                    dZ[idx] = incr[0];
                    du[idx] = incr[1];
                    dv[idx] = incr[2];
                    dw[idx] = incr[3];
                }

            // Update preconditioning of residual
            if(strcmp(itPrecon, "NO") != 0)
                applyPreconditioner(imgRows, imgCols, limit, resPCG, resCG, M);
            else
                for(int COL = limit; COL < imgCols-limit; ++COL)
                    for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                    {
                        int idx = COL*imgRows + ROW;
                        for(int i = 0; i < 4; ++i)
                            resPCG[idx][i] = resCG[idx][i];
                    }

            // beta update from all already updated residuals.
            double numBeta = 0.0;
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    numBeta += resCG[idx][0]*resPCG[idx][0] + resCG[idx][1]*resPCG[idx][1] + resCG[idx][2]*resPCG[idx][2] + resCG[idx][3]*resPCG[idx][3];
                }
            double beta = numBeta / denBeta;

            // Update direction (Preconditioned residual)
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    CG_dir(resPCG[idx], dirCG[idx], beta);
                }
        }
        // ==========================
        // Conjugate Residual Method:
        // ==========================
        else if(strcmp(itSolver, "CR") == 0)
        {
            // Updates multiplication A*p per pixel 4 unknowns
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;

                    for(int i = 0; i < 4; ++i)
                    {
                        // UP
                        stencil[i] = dirCG[COL*imgRows + (ROW-1)][i];
                        // LEFT
                        stencil[4 + i] = dirCG[(COL-1)*imgRows + ROW][i];
                        // CENTER
                        stencil[8 + i] = dirCG[COL*imgRows + ROW][i];
                        // RIGHT
                        stencil[12 + i] = dirCG[(COL+1)*imgRows + ROW][i];
                        // BOTTOM
                        stencil[16 + i] = dirCG[COL*imgRows + (ROW+1)][i];
                    }
                    mulAp(Adir[idx], aStencil[idx], stencil);
                }

            // Update preconditioning for multiplication A*p.
            if(strcmp(itPrecon, "NO") != 0)
                applyPreconditioner(imgRows, imgCols, limit, preAdir, Adir, M);
            else
                preAdir = Adir;

            // Calculate alpha value
            double numAlpha = 0.0;
            double denAlpha = 0.0;
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    for(int i = 0; i < 4; ++i)
                    {
                        // UP
                        stencil[i] = resCG[COL*imgRows + (ROW-1)][i];
                        // LEFT
                        stencil[4 + i] = resCG[(COL-1)*imgRows + ROW][i];
                        // CENTER
                        stencil[8 + i] = resCG[COL*imgRows + ROW][i];
                        // RIGHT
                        stencil[12 + i] = resCG[(COL+1)*imgRows + ROW][i];
                        // BOTTOM
                        stencil[16 + i] = resCG[COL*imgRows + (ROW+1)][i];
                    }
                    mulAp(Ares[idx], aStencil[idx], stencil);
                    numAlpha += resCG[idx][0]*Ares[idx][0] + resCG[idx][1]*Ares[idx][1] + resCG[idx][2]*Ares[idx][2] + resCG[idx][3]*Ares[idx][3];
                    denAlpha += Adir[idx][0]*preAdir[idx][0] + Adir[idx][1]*preAdir[idx][1] + Adir[idx][2]*preAdir[idx][2] + Adir[idx][3]*preAdir[idx][3];
                }
            double denBeta = numAlpha;
            double alpha = numAlpha / denAlpha;

            // Calculation of updated increments and residual
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    incr[0] = dZ[idx];
                    incr[1] = du[idx];
                    incr[2] = dv[idx];
                    incr[3] = dw[idx];

                    CG_incr_res(incr, resCG[idx], dirCG[idx], preAdir[idx], alpha);
                    for(int i = 0; i < 4; ++i)
                        residual[idx][i] = resCG[idx][i];

                    dZ[idx] = incr[0];
                    du[idx] = incr[1];
                    dv[idx] = incr[2];
                    dw[idx] = incr[3];
                }

            // beta update from all updated residuals
            double numBeta = 0.0;
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    for(int i = 0; i < 4; ++i)
                    {
                        // UP
                        stencil[i] = resCG[COL*imgRows + (ROW-1)][i];
                        // LEFT
                        stencil[4 + i] = resCG[(COL-1)*imgRows + ROW][i];
                        // CENTER
                        stencil[8 + i] = resCG[COL*imgRows + ROW][i];
                        // RIGHT
                        stencil[12 + i] = resCG[(COL+1)*imgRows + ROW][i];
                        // BOTTOM
                        stencil[16 + i] = resCG[COL*imgRows + (ROW+1)][i];
                    }
                    mulAp(Ares[idx], aStencil[idx], stencil);
                    numBeta += resCG[idx][0]*Ares[idx][0] + resCG[idx][1]*Ares[idx][1] + resCG[idx][2]*Ares[idx][2] + resCG[idx][3]*Ares[idx][3];
                }
            double beta = numBeta / denBeta;

            // Update direction (Not preconditioned residual)
            for(int COL = limit; COL < imgCols-limit; ++COL)
                for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                {
                    int idx = COL*imgRows + ROW;
                    CG_dir(resCG[idx], dirCG[idx], beta);
                }


            // Calculate Residual (preconditioned is very very small)
            if(strcmp(itPrecon, "NO") != 0)
                for(int COL = limit; COL < imgCols-limit; ++COL)
                    for(int ROW = limit; ROW < imgRows-limit; ++ROW)
                    {
                        int idx = COL*imgRows + ROW;
                        int idxUp = COL*imgRows + ROW-1;
                        int idxLeft = (COL-1)*imgRows + ROW;
                        int idxRight = (COL+1)*imgRows + ROW;
                        int idxDown = COL*imgRows + ROW+1;

                        stencil[0] = dZ[idxUp]; stencil[1] = du[idxUp]; stencil[2] = dv[idxUp]; stencil[3] = dw[idxUp];
                        stencil[4] = dZ[idxLeft]; stencil[5] = du[idxLeft]; stencil[6] = dv[idxLeft]; stencil[7] = dw[idxLeft];
                        stencil[8] = dZ[idx]; stencil[9] = du[idx]; stencil[10] = dv[idx]; stencil[11] = dw[idx];
                        stencil[12] = dZ[idxRight]; stencil[13] = du[idxRight]; stencil[14] = dv[idxRight]; stencil[15] = dw[idxRight];
                        stencil[16] = dZ[idxDown]; stencil[17] = du[idxDown]; stencil[18] = dv[idxDown]; stencil[19] = dw[idxDown];

                        mulAp(residual[idx], aStencil[idx], stencil);
                        residual[idx][0] = b[idx + rowPxl*0] - residual[idx][0];
                        residual[idx][1] = b[idx + rowPxl*1] - residual[idx][1];
                        residual[idx][2] = b[idx + rowPxl*2] - residual[idx][2];
                        residual[idx][3] = b[idx + rowPxl*3] - residual[idx][3];
                    }

        }
        // ==========================
        // Successive over-relaxation:
        // ==========================
        else if(strcmp(itSolver, "SOR") == 0)
        {
            // 4 subiterations from top-bottom, left-right and viceversa
            int initRow, endRow , stepRow, initCol, endCol, stepCol;
            if(scanMode == 0)
            {
                initRow = 0; endRow = imgRows; stepRow = 1;
                initCol = 0; endCol = imgCols; stepCol = 1;
            }
            else if(scanMode == 1)
            {
                initRow = imgRows-1; endRow = -1; stepRow = -1;
                initCol = imgCols-1; endCol = -1; stepCol = -1;
            }
            else if(scanMode == 2)
            {
                initRow = 0; endRow = imgRows; stepRow = 1;
                initCol = imgCols-1; endCol = -1; stepCol = -1;
            }
            else if (scanMode == 3)
            {
                initRow = imgRows-1; endRow = -1; stepRow = -1;
                initCol = 0; endCol = imgCols; stepCol = 1;
            }

            for(int col = initCol; col != endCol; col += stepCol)
                for(int row = initRow; row != endRow; row += stepRow)
                {
                    int idx = col*imgRows + row;

                    // Special treatment for borders: Mirroring - Neumann boundary conditions.
                    // (take pixel value of adjacent pxl's other side)
                    if(row <= limit || row >= imgRows-limit-1 || col <= limit || col >= imgCols-limit-1)
                    {
                        Increments border = treatBorders(imgRows, imgCols, row, col, dZ, du, dv, dw, limit+1);
                        dZ[idx] = border.dZ;
                        du[idx] = border.du;
                        dv[idx] = border.dv;
                        dw[idx] = border.dw;
                        continue;
                    }

                    double* coef = aStencil[idx];
                    // div_dZ
                    div_d[0] = coef[0*4+0]*dZ[col*imgRows + row-1] + coef[4*4+0]*dZ[(col*imgRows-imgRows) + row] +
                            coef[12*4+0]*dZ[(col*imgRows+imgRows) + row] + coef[16*4+0]*dZ[col*imgRows + row+1];
                    // div_du
                    div_d[1] = coef[1*4+1]*du[col*imgRows + row-1] + coef[5*4+1]*du[(col*imgRows-imgRows) + row] +
                            coef[13*4+1]*du[(col*imgRows+imgRows) + row] + coef[17*4+1]*du[col*imgRows + row+1];
                    // div_dv
                    div_d[2] = coef[2*4+2]*dv[col*imgRows + row-1] + coef[6*4+2]*dv[(col*imgRows-imgRows) + row] +
                            coef[14*4+2]*dv[(col*imgRows+imgRows) + row] + coef[18*4+2]*dv[col*imgRows + row+1];
                    // div_dw
                    div_d[3] = coef[3*4+3]*dw[col*imgRows + row-1] + coef[7*4+3]*dw[(col*imgRows-imgRows) + row] +
                            coef[15*4+3]*dw[(col*imgRows+imgRows) + row] + coef[19*4+3]*dw[col*imgRows + row+1];

                    incr[0] = dZ[idx];
                    incr[1] = du[idx];
                    incr[2] = dv[idx];
                    incr[3] = dw[idx];

                    Increments pxlSolver = SOR(imgRows, imgCols, a, b, row, col, incr, div_d, wSOR);

                    // Residual of current pixel.
                    for(int i = 0; i < 4; ++i)
                        bResidual[i] = b[idx + rowPxl*i] - div_d[i];
                    residualSystem(imgRows, imgCols, residual[idx], a, row, col, bResidual, incr);

                    dZ[idx] = pxlSolver.dZ;
                    du[idx] = pxlSolver.du;
                    dv[idx] = pxlSolver.dv;
                    dw[idx] = pxlSolver.dw;
                }

            scanMode = (scanMode + 1) % 4;
        }

        // Check convergence with previous and current iteration residuals.
        int limitBorder = limit;
        if(strcmp(itSolver, "SOR") == 0)
            limitBorder = limit+1;
        double evol = 0.0;
        for(int COL = limitBorder; COL < imgCols-limitBorder; ++COL)
            for(int ROW = limitBorder; ROW < imgRows-limitBorder; ++ROW)
                for(int incr = 0; incr < 4; ++incr)
                    evol += residual[imgRows*COL + ROW][incr]*residual[imgRows*COL + ROW][incr];
        double relativeResidual = sqrt(evol);
        relativeResidual = (relativeResidual / valueB)*100.0;
        evolResidual[it] = relativeResidual;
        convergence = (relativeResidual < errorSolver);
        it++;
    }

    // mexPrintf("END\n");

    // Special treatment for borders: Neumann bound. cond.
    // (take pixel value of adjacent pxl's other side)
    if(strcmp(itSolver, "SOR") != 0)
    {
        // Top & Bottom
        for(int COL = limit; COL < imgCols-limit; ++COL)
        {
            for(int ROW = 0; ROW < limit; ++ROW)
            {
                dZ[imgRows*COL + ROW] = dZ[imgRows*COL + limit];
                du[imgRows*COL + ROW] = du[imgRows*COL + limit];
                dv[imgRows*COL + ROW] = dv[imgRows*COL + limit];
                dw[imgRows*COL + ROW] = dw[imgRows*COL + limit];
            }
            for(int ROW = imgRows-limit; ROW < imgRows; ++ROW)
            {
                dZ[imgRows*COL + ROW] = dZ[imgRows*COL + imgRows-limit-1];
                du[imgRows*COL + ROW] = du[imgRows*COL + imgRows-limit-1];
                dv[imgRows*COL + ROW] = dv[imgRows*COL + imgRows-limit-1];
                dw[imgRows*COL + ROW] = dw[imgRows*COL + imgRows-limit-1];
            }
        }
        // Left & Right
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            for(int COL = 0; COL < limit; ++COL)
            {
                dZ[imgRows*COL + ROW] = dZ[imgRows*limit + ROW];
                du[imgRows*COL + ROW] = du[imgRows*limit + ROW];
                dv[imgRows*COL + ROW] = dv[imgRows*limit + ROW];
                dw[imgRows*COL + ROW] = dw[imgRows*limit + ROW];
            }
            for(int COL = imgCols-limit; COL < imgCols; ++COL)
            {
                dZ[imgRows*COL + ROW] = dZ[imgRows*(imgCols-limit-1) + ROW];
                du[imgRows*COL + ROW] = du[imgRows*(imgCols-limit-1) + ROW];
                dv[imgRows*COL + ROW] = dv[imgRows*(imgCols-limit-1) + ROW];
                dw[imgRows*COL + ROW] = dw[imgRows*(imgCols-limit-1) + ROW];
            }
        }
    }
    // For SOR make the mean of 4-iteration group.
    else
    {
        for(int i = 0; i < it; i+=4)
        {
            if(i+3 < it)
            {
                evolResidual[i] = evolResidual[i+3]; evolResidual[i+1] = evolResidual[i+3]; evolResidual[i+2] = evolResidual[i+3];
            }
            else if(i+1 >= it)
            {
                break;
            }
            else if(i+2 >= it)
            {
                evolResidual[i] = evolResidual[i+1];
            }
            else
            {
                evolResidual[i] = evolResidual[i+2];
                evolResidual[i+1] = evolResidual[i+2];
            }
        }
    }

    // Release memory

    for(int COL = limit; COL < imgCols-limit; ++COL)
        for(int ROW = limit; ROW < imgRows-limit; ++ROW)
        {
            delete [] aStencil[imgRows*COL + ROW];
            delete [] residual[imgRows*COL + ROW];
        }

    delete [] aStencil;
    delete [] residual;

    // CG & CR
    if(strcmp(itSolver, "CG") == 0 || strcmp(itSolver, "CR") == 0)
    {
        for(int COL = 0; COL < imgCols; ++COL)
            for(int ROW = 0; ROW < imgRows; ++ROW)
            {
                delete [] resCG[imgRows*COL + ROW];
                delete [] dirCG[imgRows*COL + ROW];
                delete [] resPCG[imgRows*COL + ROW];
                delete [] Adir[imgRows*COL + ROW];
                delete [] Ares[imgRows*COL + ROW];
            }
        delete [] resCG;
        delete [] dirCG;
        delete [] resPCG;
        delete [] Adir;
        delete [] Ares;

        if(strcmp(itPrecon, "NO") != 0)
        {
            for(int COL = 0; COL < imgCols; ++COL)
                for(int ROW = 0; ROW < imgRows; ++ROW)
                {
                    int idx = imgRows*COL + ROW;
                    delete [] M[idx];
                }

            delete [] M;

            if(strcmp(itSolver, "CR") == 0)
            {
                for(int COL = 0; COL < imgCols; ++COL)
                    for(int ROW = 0; ROW < imgRows; ++ROW)
                    {
                        int idx = imgRows*COL + ROW;
                        delete [] preAdir[idx];
                    }
                    delete [] preAdir;
            }
        }
    }
}

// Solves a system of linear equations with SOR algorithm.
// INPUT PARAMETERS:
// - 4x4 A matrix with unknown dep}ent konstants
// - 1x4 b matrix with unknown indep}ent konstants
Increments SOR(int imgRows, int imgCols, double* a, double* b, int rowMat, int colMat, double* incr, double* div_d, double w)
{
    int idx = colMat*imgRows + rowMat;
    int rowPxl = imgRows*imgCols;
    int colPxl = imgRows*imgCols*4;

    // Storage of previous value to measure level of convergence.
    for(int row = 0; row < 4; ++row)
    {
        // Reset sum next row.
        double sigma = 0.0;
        for(int col = 0; col < 4; ++col)
        {
            if(row != col)
                sigma = sigma + a[idx + rowPxl*row + colPxl*col] * incr[col];
        }
        incr[row] = (1 - w)*incr[row] + w * ( (b[idx + rowPxl*row] - sigma - div_d[row]) / a[idx + rowPxl*row + colPxl*row] );
    }

    Increments incrOut;
    incrOut.dZ = incr[0];
    incrOut.du = incr[1];
    incrOut.dv = incr[2];
    incrOut.dw = incr[3];
    return incrOut;
}

// Solves a system of linear equations with Conjugate Gradient method.
// INPUT PARAMETERS:
// - 4x4 A matrix with unknown dependent konstants
// - 1x4 b matrix with unknown independent konstants
void CG_incr_res(double* incr, double* r, double* p, double* Ap, double alpha)
{
    // Update of increments
    incr[0] += alpha*p[0];
    incr[1] += alpha*p[1];
    incr[2] += alpha*p[2];
    incr[3] += alpha*p[3];

    // Update of the residual
    r[0] -= alpha*Ap[0];
    r[1] -= alpha*Ap[1];
    r[2] -= alpha*Ap[2];
    r[3] -= alpha*Ap[3];
}

void CG_dir(double* r, double* d, double beta)
{
    d[0] = r[0] + beta*d[0];
    d[1] = r[1] + beta*d[1];
    d[2] = r[2] + beta*d[2];
    d[3] = r[3] + beta*d[3];
}

// Private function: Residual of approximation
// INPUT PARAMETERS:
// - 4x4 A matrix with unknown dep}ent konstants
// - 1x4 b matrix with unknown indep}ent konstants
// - 1x4 x matrix with approximated unknowns
void residualSystem(int imgRows, int imgCols, double* residual, double* a, int row, int col, double* b, double* x)
{
    int idx = col*imgRows + row;
    int rowPxl = imgRows*imgCols;
    int colPxl = imgRows*imgCols*4;

    // result = norm(b - a*x)
    residual[0] = b[0] - (a[idx + rowPxl*0 + colPxl*0]*x[0] + a[idx + rowPxl*0 + colPxl*1]*x[1] + a[idx + rowPxl*0 + colPxl*2]*x[2] + a[idx + rowPxl*0 + colPxl*3]*x[3]);
    residual[1] = b[1] - (a[idx + rowPxl*1 + colPxl*0]*x[0] + a[idx + rowPxl*1 + colPxl*1]*x[1] + a[idx + rowPxl*1 + colPxl*2]*x[2] + a[idx + rowPxl*1 + colPxl*3]*x[3]);
    residual[2] = b[2] - (a[idx + rowPxl*2 + colPxl*0]*x[0] + a[idx + rowPxl*2 + colPxl*1]*x[1] + a[idx + rowPxl*2 + colPxl*2]*x[2] + a[idx + rowPxl*2 + colPxl*3]*x[3]);
    residual[3] = b[3] - (a[idx + rowPxl*3 + colPxl*0]*x[0] + a[idx + rowPxl*3 + colPxl*1]*x[1] + a[idx + rowPxl*3 + colPxl*2]*x[2] + a[idx + rowPxl*3 + colPxl*3]*x[3]);
}

// Private function: Assign borders incr's value from adjacent pixels.
// INPUT PARAMETERS:
// - Row in the img matrix
// - Col in the img matrix
// - matrix of increments struct
Increments treatBorders(int imgRows, int imgCols, int ROW, int COL, double* dZ, double* du, double* dv, double* dw, int limit)
{
    Increments incr;
    incr.initZero();

    if(ROW < limit)
    {
        // Top left corner.
        if(COL < limit)
        {
            incr.dZ = dZ[limit*imgRows + limit];
            incr.du = du[limit*imgRows + limit];
            incr.dv = dv[limit*imgRows + limit];
            incr.dw = dw[limit*imgRows + limit];
        }
        // Top right corner.
        else if(COL >= imgCols-limit)
        {
            incr.dZ = dZ[(imgCols-limit-1)*imgRows + limit];
            incr.du = du[(imgCols-limit-1)*imgRows + limit];
            incr.dv = dv[(imgCols-limit-1)*imgRows + limit];
            incr.dw = dw[(imgCols-limit-1)*imgRows + limit];
        }
        // Top border (no corner).
        else
        {
            incr.dZ = dZ[COL*imgRows + limit];
            incr.du = du[COL*imgRows + limit];
            incr.dv = dv[COL*imgRows + limit];
            incr.dw = dw[COL*imgRows + limit];
        }
    }
    else if(ROW >= imgRows-limit)
    {
        // Bottom left corner.
        if(COL < limit)
        {
            incr.dZ = dZ[limit*imgRows + imgRows-limit-1];
            incr.du = du[limit*imgRows + imgRows-limit-1];
            incr.dv = dv[limit*imgRows + imgRows-limit-1];
            incr.dw = dw[limit*imgRows + imgRows-limit-1];
        }
        // Bottom right corner.
        else if(COL >= imgCols-limit)
        {
            incr.dZ = dZ[(imgCols-limit-1)*imgRows + imgRows-limit-1];
            incr.du = du[(imgCols-limit-1)*imgRows + imgRows-limit-1];
            incr.dv = dv[(imgCols-limit-1)*imgRows + imgRows-limit-1];
            incr.dw = dw[(imgCols-limit-1)*imgRows + imgRows-limit-1];
        }
        // Bottom border (no corner).
        else
        {
            incr.dZ = dZ[COL*imgRows + imgRows-limit-1];
            incr.du = du[COL*imgRows + imgRows-limit-1];
            incr.dv = dv[COL*imgRows + imgRows-limit-1];
            incr.dw = dw[COL*imgRows + imgRows-limit-1];
        }
    }
    else
    {
        // Left side (no corner).
        if(COL < limit)
        {
            incr.dZ = dZ[limit*imgRows + ROW];
            incr.du = du[limit*imgRows + ROW];
            incr.dv = dv[limit*imgRows + ROW];
            incr.dw = dw[limit*imgRows + ROW];
        }
        // Right side (no corner).
        else if(COL >= imgCols-limit)
        {
            incr.dZ = dZ[(imgCols-limit-1)*imgRows + ROW];
            incr.du = du[(imgCols-limit-1)*imgRows + ROW];
            incr.dv = dv[(imgCols-limit-1)*imgRows + ROW];
            incr.dw = dw[(imgCols-limit-1)*imgRows + ROW];
        }
    }

    return incr;
}

void mulAp(double* Ap, double* A, double* p)
{
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 1; ++j)
        {
            Ap[j*4 + i] = 0.0;
            for(int k = 0; k < 20; ++k)
            {
                Ap[j*4 + i] += A[k*4 + i]*p[j*20 + k];
            }
        }
}

/* Preconditioning */
void applyPreconditioner(int imgRows, int imgCols, int limit, double** z, double** r, double** M)
{
    // Go from (Mz = r) to (z = M^-1 r)
    // M = (D+L)D^-1(D+U) -> M = (D+L)(I + D^-1 U)
    // Hence: M = LU -> (1) L = (D+L) (2) U = (I + D^-1 U)

    // (1) Loop lower triangular part of the matrix (forward sweep)
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = imgRows*COL + ROW;
            z[idx][0] = 0.0; z[idx][1] = 0.0; z[idx][2] = 0.0; z[idx][3] = 0.0;

            // Treat boundaries: all will be 0.
            if(ROW < limit || ROW >= imgRows-limit || COL < limit || COL >= imgCols-limit)
                continue;

            for(int incr = 0; incr < 4; ++incr)
            {
                double nonZeros = 0.0;

                // Top Pxl increments
                nonZeros += M[idx][incr*4 + incr]*z[imgRows*COL + ROW-1][incr];
                // Left Pxl increments
                nonZeros += M[idx][(incr+4)*4 + incr]*z[imgRows*(COL-1) + ROW][incr];
                // Center pixel increments
                for(int pxl = 0; pxl < incr; ++pxl)
                    nonZeros += M[idx][(8+pxl)*4 + incr]*z[idx][pxl];

                // Isolate and solve half-preconditioned z
                z[idx][incr] = (r[idx][incr] - nonZeros) / M[idx][(8+incr)*4 + incr];
            }
        }

    // (2) Loop lower triangular part of the matrix (forward sweep)
    // Not necessary to loop the borders anymore.
    for(int COL = imgCols-limit-1; COL >= limit ; --COL)
        for(int ROW = imgRows-limit-1; ROW >= limit ; --ROW)
        {
            int idx = COL*imgRows + ROW;

            for(int incr = 3; incr >= 0; --incr)
            {
                double nonZeros = 0.0;

                // Right Pxl increments
                nonZeros += M[idx][(incr+12)*4 + incr]*z[imgRows*(COL+1) + ROW][incr];
                // Bottom Pxl increments
                nonZeros += M[idx][(incr+16)*4 + incr]*z[imgRows*COL + ROW+1][incr];
                // Center increments
                for(int pxl = 3; pxl > incr; --pxl)
                    nonZeros += M[idx][(8+pxl)*4 + incr]*z[idx][pxl];

                // Isolate and solve full-preconditioned z
                z[idx][incr] -= nonZeros / M[idx][(8+incr)*4 + incr];
            }
        }
}

void computePreconditioner(double** M, int imgRows, int imgCols, int limit, char* preconditioner, double** A)
{
    // We start preconditioning storing whole A efficiently.
    // This is necessary for all kinds of preconditions.
    // SSOR won't need any further step.
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            if(ROW < limit || ROW >= imgRows-limit || COL < limit || COL >= imgCols-limit)
            {
                for(int coef = 0; coef < 4*20; ++coef)
                    M[imgRows*COL + ROW][coef] = 0.0;
                continue;
            }
            else
                for(int coef = 0; coef < 4*20; ++coef)
                    M[imgRows*COL + ROW][coef] = A[imgRows*COL + ROW][coef];
        }

    // Calculate preconditioner for incomplete factorization cases.
    if(strcmp(preconditioner, "ILU") == 0)
        computeILU(imgRows, imgCols, limit, M);
}

void computeILU(int imgRows, int imgCols, int limit, double** A)
{
    // ILU factorization, IKJ version.
    // For all rows...
    double pivot = 0.0;
    for(int COL = limit; COL < imgCols-limit; ++COL)
        for(int ROW = limit; ROW < imgRows-limit; ++ROW)
        {
            int idx = imgRows*COL + ROW;
            for(int incr = 0; incr < 4; ++incr)
            {
                // - Left pixel:
                // Calculate pivot = A(i,k) / A(k,k);
                if(A[idx][(incr+4)*4 + incr] != 0)
                {
                    pivot = A[idx][(incr+4)*4 + incr] / A[imgRows*(COL-1) + ROW][(8+incr)*4 + incr];
                    // Apply to top pixel
                    A[idx][(incr+8)*4 + incr] -= pivot * A[imgRows*(COL-1) + ROW][(12+incr)*4 + incr];
                }

                // - Top pixel:
                // Calculate pivot = A(i,k) / A(k,k);
                if(A[idx][incr*4 + incr] != 0)
                {
                    pivot = A[idx][incr*4 + incr] / A[imgRows*COL + ROW-1][(8+incr)*4 + incr];
                    // Apply to center pixels
                    A[idx][(incr+8)*4 + incr] -= pivot * A[imgRows*(COL) + ROW-1][(16+incr)*4 + incr];
                }

                // - Center pixels up to current 'incr'
                for(int pivotIncr = 0; pivotIncr < incr; ++pivotIncr)
                {
                    // Calculate pivot = A(i,k) / A(k,k);
                    if(A[idx][(pivotIncr+8)*4 + incr] != 0)
                    {
                        pivot = A[idx][(pivotIncr+8)*4 + incr] / A[imgRows*COL + ROW][(8+pivotIncr)*4 + pivotIncr];

                        // Apply to center pixels
                        for(int j = pivotIncr+1; j < 4; ++j)
                            A[idx][(j+8)*4 + incr] -= pivot * A[imgRows*COL + ROW][(8+j)*4 + pivotIncr];
                    }
                }
            }
        }
}


