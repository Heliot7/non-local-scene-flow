#include <math.h>
#include <algorithm>
#include <limits>
#include <cstring>

#include "mex.h" // for: mexPrintf(...)

#include "libMath/matrix.hpp"
#include "differentiation.hpp"
#include "disparity.hpp"
#include "colour.hpp"

using namespace std;

// Estimates an intial disparity with segmentation as a prior and voting.
void computeLocalDisparities(int* dispL, double* scoresL, int* dispR, double* scoresR, int imgRows, int imgCols, int* segL, int numSegL, int* segR, int numSegR,
                             unsigned char* imageL, unsigned char* imageR, int imgChannels, int maxDisparity, double dispLambda, int baseline, char* mCost)
{
    // Initialise matching costs for both images.
    double* matchCostL  = new double[imgRows*imgCols*maxDisparity];
    double* matchCostR  = new double[imgRows*imgCols*maxDisparity];
    for(int idx = 0; idx < imgRows*imgCols*maxDisparity; ++idx)
    {
        matchCostL[idx] = 0.0;
        matchCostR[idx] = 0.0;
    }

    if(strcmp(mCost, "SAD") == 0)
    {
        // Intensity/colour channel costs.
        for(int CHANNEL = 0; CHANNEL < imgChannels; ++CHANNEL)
            calculateIndividualCosts(imgRows, imgCols, maxDisparity, CHANNEL, imageL, matchCostL, imageR, matchCostR);

        // Calculate gradient costs (X and Y dimensions are channels 1 and 2, then).
        double* imgLxy = new double[imgRows*imgCols*2];
        double* imgRxy = new double[imgRows*imgCols*2];

        gradient(imgRows, imgCols, imageL, imgLxy);
        gradient(imgRows, imgCols, imageR, imgRxy);

        // TO DO!
        calculateIndividualCosts(imgRows, imgCols, maxDisparity, 0, imgLxy, matchCostL, imgRxy, matchCostR); // X
        calculateIndividualCosts(imgRows, imgCols, maxDisparity, 1, imgLxy, matchCostL, imgRxy, matchCostR); // Y

        delete [] imgLxy; delete [] imgRxy;
    }
    else if(strcmp(mCost, "CENSUS") == 0)
    {
        unsigned char** censusL, **censusR;
        censusL = new unsigned char*[imgRows*imgCols];
        censusR = new unsigned char*[imgRows*imgCols];
        for(int idx = 0; idx < imgRows*imgCols; ++idx)
        {
            censusL[idx] = new unsigned char[9*7];
            censusR[idx] = new unsigned char[9*7];
            for(int c = 0; c < 9*7; ++c)
            {
                censusL[idx][c] = 0;
                censusR[idx][c] = 0;
            }
        }

        unsigned char *imgGreyL = new unsigned char[imgRows*imgCols];
        unsigned char *imgGreyR = new unsigned char[imgRows*imgCols];
        if(imgChannels == 3)
        {
            rgbToGreyscale(imgRows, imgCols, imgGreyL, imageL);
            rgbToGreyscale(imgRows, imgCols, imgGreyR, imageR);
        }
        else if(imgChannels == 1)
        {
            imgGreyL = imageL;
            imgGreyR = imageR;
        }

        imageCensus(censusL, imgRows, imgCols, imgGreyL);
        imageCensus(censusR, imgRows, imgCols, imgGreyR);
        censusCost(imgRows, imgCols, maxDisparity, censusL, matchCostL, censusR, matchCostR);

        for(int idx = 0; idx < imgRows*imgCols; ++idx)
        {
            delete [] censusL[idx];
            delete [] censusR[idx];
        }
        delete [] censusL;
        delete [] censusR;
    }
    else // BILSUB
    {

    }

    // Window size aggregation
    int w = (int)(0.10 * imgCols);
    // mexPrintf("Window size aggregation: %d\n", w);

    // Create auxiliar score for images 1 and 2
    for(int i = 0; i < imgRows*imgCols*maxDisparity; ++i)
    {
        scoresL[i] = numeric_limits<double>::infinity();
        scoresR[i] = numeric_limits<double>::infinity();
    }
    // Compute aggregated costs
    double* costPxlSegmentL = new double[imgRows*imgCols];
    double* costPxlSegmentR = new double[imgRows*imgCols];
    double* costPxlWindowL = new double[imgRows*imgCols];
    double* costPxlWindowR = new double[imgRows*imgCols];
    double costSegmentL[numSegL];
    double costSegmentR[numSegR];
    for(int DISP = 0; DISP < maxDisparity; ++DISP)
    {
        for(int i = 0; i < numSegL; ++i)
            costSegmentL[i] = 0.0;
        for(int i = 0; i < numSegR; ++i)
            costSegmentR[i] = 0.0;

        for(int i = 0; i < imgRows*imgCols; ++i)
        {
            costPxlSegmentL[i] = 0.0;
            costPxlSegmentR[i] = 0.0;
            costPxlWindowL[i] = 0.0;
            costPxlWindowR[i] = 0.0;
        }

        // 1st Sweep: Through all rows...
        for(int ROW = 0; ROW < imgRows; ++ROW)
            for(int COL = (int)-floor(w/2); COL < imgCols+(int)floor(w/2); ++COL)
            {
                double borderL = 0.0;
                double borderR = 0.0;

                // Check pxl borders to do or not do.
                int rightWinPxl = COL + (int)floor(w/2);
                if(rightWinPxl < imgCols)
                {
                    borderL += matchCostL[imgRows*rightWinPxl + ROW + imgRows*imgCols*DISP];
                    costSegmentL[segL[imgRows*rightWinPxl + ROW]-1] += matchCostL[imgRows*rightWinPxl + ROW + imgRows*imgCols*DISP];

                    borderR += matchCostR[imgRows*rightWinPxl + ROW + imgRows*imgCols*DISP];
                    costSegmentR[segR[imgRows*rightWinPxl + ROW]-1] += matchCostR[imgRows*rightWinPxl + ROW + imgRows*imgCols*DISP];
                }

                int leftWinPxl = COL - (int)floor(w/2);
                if(leftWinPxl >= 0)
                {
                    borderL -= matchCostL[imgRows*leftWinPxl + ROW + imgRows*imgCols*DISP];
                    costSegmentL[segL[imgRows*leftWinPxl + ROW]-1] -= matchCostL[imgRows*leftWinPxl + ROW + imgRows*imgCols*DISP];

                    borderR -= matchCostR[imgRows*leftWinPxl + ROW + imgRows*imgCols*DISP];
                    costSegmentR[segR[imgRows*leftWinPxl + ROW]-1] -= matchCostR[imgRows*leftWinPxl + ROW + imgRows*imgCols*DISP];
                }

                // Assign the aggregated cost
                if(COL >= 0 && COL < imgCols)
                {
                    costPxlSegmentL[imgRows*COL + ROW] = costSegmentL[segL[imgRows*COL + ROW]-1];
                    costPxlSegmentR[imgRows*COL + ROW] = costSegmentR[segR[imgRows*COL + ROW]-1];
                }

                if(COL <= 0)
                {
                    costPxlWindowL[ROW] += borderL;
                    costPxlWindowR[ROW] += borderR;
                }
                else if(COL > 0 && COL < imgCols)
                {
                    costPxlWindowL[imgRows*COL + ROW] = costPxlWindowL[imgRows*(COL-1) + ROW] + borderL;
                    costPxlWindowR[imgRows*COL + ROW] = costPxlWindowR[imgRows*(COL-1) + ROW] + borderR;
                }
            }

        // 2nd Sweep: Through all cols...
        for(int i = 0; i < numSegL; ++i)
            costSegmentL[i] = 0.0;
        for(int i = 0; i < numSegR; ++i)
            costSegmentR[i] = 0.0;
        double tL = 0.0;
        double tR = 0.0;
        for(int COL = 0; COL < imgCols; ++COL)
            for(int ROW = (int)-floor(w/2); ROW < imgRows+(int)floor(w/2); ++ROW)
            {
                int downWinPxl = ROW + (int)floor(w/2);
                int upWinPxl= ROW - (int)floor(w/2);

                // Segment info
                if(downWinPxl < imgRows)
                {
                    costSegmentL[segL[imgRows*COL + downWinPxl]-1] += costPxlSegmentL[imgRows*COL + downWinPxl];
                    costSegmentR[segR[imgRows*COL + downWinPxl]-1] += costPxlSegmentR[imgRows*COL + downWinPxl];
                }

                if(upWinPxl >= 0)
                {
                    costSegmentL[segL[imgRows*COL + upWinPxl]-1] -= costPxlSegmentL[imgRows*COL + upWinPxl];
                    costSegmentR[segR[imgRows*COL + upWinPxl]-1] -= costPxlSegmentR[imgRows*COL + upWinPxl];
                }

                // Window info
                if(downWinPxl < imgRows)
                {
                    tL += costPxlWindowL[imgRows*COL + downWinPxl];
                    tR += costPxlWindowR[imgRows*COL + downWinPxl];
                }
                if(upWinPxl >= 0)
                {
                    tL -= costPxlWindowL[imgRows*COL + upWinPxl];
                    tR -= costPxlWindowR[imgRows*COL + upWinPxl];
                }

                // Append scores
                if(ROW >= 0 && ROW < imgRows)
                {
                    scoresL[imgRows*COL + ROW + imgRows*imgCols*DISP] = dispLambda*(tL - costSegmentL[segL[imgRows*COL + ROW]-1]) + costSegmentL[segL[imgRows*COL + ROW]-1];
                    scoresR[imgRows*COL + ROW + imgRows*imgCols*DISP] = dispLambda*(tR - costSegmentR[segR[imgRows*COL + ROW]-1]) + costSegmentR[segR[imgRows*COL + ROW]-1];
                }
            }
    }

    // Assign new score to below baseline pixels in those disparities out of bounds
    for(int COL = 0; COL < baseline; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
            for(int DISP = COL; DISP < maxDisparity; ++DISP)
            {
                scoresL[imgRows*COL + ROW + imgRows*imgCols*DISP] = scoresL[imgRows*baseline + ROW + imgRows*imgCols*DISP];
                scoresR[imgRows*(imgCols-1-COL) + ROW + imgRows*imgCols*DISP] = scoresR[imgRows*(imgCols-1-baseline) + ROW + imgRows*imgCols*DISP];
            }

    // Choose initial best matching scores (no plane dependance).
    double bestScore, newScore;
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            bestScore = scoresL[imgRows*COL + ROW];
            int bestDisp = 0;
            for(int DISP = 1; DISP < maxDisparity; ++DISP)
            {
                newScore = scoresL[imgRows*COL + ROW + imgRows*imgCols*DISP];
                if(newScore < bestScore)
                {
                    bestScore = newScore;
                    bestDisp = DISP;
                }
            }
            dispL[imgRows*COL + ROW] = (double)bestDisp + 1.0;
        }
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            bestScore = scoresR[imgRows*COL + ROW];
            int bestDisp = 0;
            for(int DISP = 1; DISP < maxDisparity; ++DISP)
            {
                newScore = scoresR[imgRows*COL + ROW + imgRows*imgCols*DISP];
                if(newScore < bestScore)
                {
                    bestScore = newScore;
                    bestDisp = DISP;
                }
            }
            dispR[imgRows*COL + ROW] = (double)bestDisp + 1.0;
        }

    delete [] matchCostL; delete [] matchCostR;
    delete [] costPxlSegmentL; delete [] costPxlSegmentR;
    delete [] costPxlWindowL; delete [] costPxlWindowR;
}

// Private function: Set individual matching costs among all pixel candidates and their disparities.
template <class T>
void calculateIndividualCosts(int imgRows, int imgCols, int maxDisparity, int channel, T* imageL, double* costL, T* imageR, double* costR)
{
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            double bigValueL = 0.0;
            double bigValueR = 0.0;
            for(int DISP = 1; DISP <= maxDisparity; ++DISP)
            {
                if(COL-DISP >= 0)
                {
                    double newValue = fabs(imageL[imgRows*COL + ROW + imgRows*imgCols*channel] - imageR[imgRows*(COL-DISP) + ROW + imgRows*imgCols*channel]);
                    if(newValue > bigValueL)
                        bigValueL = newValue;
                    costL[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += newValue;
                }
                else
                    costL[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += bigValueL;

                if(COL+DISP < imgCols)
                {
                    double newValue = fabs(imageL[imgRows*(COL+DISP) + ROW + imgRows*imgCols*channel] - imageR[imgRows*COL + ROW + imgRows*imgCols*channel]);
                    if(newValue > bigValueR)
                        bigValueR = newValue;
                    costR[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += newValue;
                }
                else
                    costR[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += bigValueR;
            }
        }
}

// Computes the census of the current image
void imageCensus(unsigned char** census, int imgRows, int imgCols, unsigned char* image)
{
    for(int col = 0; col < imgCols; ++col)
        for(int row = 0; row < imgRows; ++row)
        {
            int posCensus = 0;
            for(int colWin = col-3; colWin < col+4; ++colWin)
                for(int rowWin = row-4; rowWin < row+4; ++rowWin)
                {
                    posCensus++;

                    // Ignore comparisons out of bounds
                    if(rowWin < 0 || rowWin >= imgRows || colWin < 0 || colWin >= imgCols)
                        continue;

                    // Compare neighbours to update census filter
                    if(image[colWin*imgRows+rowWin] < image[col*imgRows+row])
                        census[col*imgRows+row][posCensus] = 1;
                }
        }
}

void censusCost(int imgRows, int imgCols, int maxDisparity, unsigned char** censusL, double* costL, unsigned char** censusR, double* costR)
{
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int bigValueL = 0;
            int bigValueR = 0;
            for(int DISP = 1; DISP <= maxDisparity; ++DISP)
            {
                if(COL-DISP >= 0)
                {
                    int hammingDist = censusHammingDistance(censusL[imgRows*COL + ROW], censusR[imgRows*(COL-DISP) + ROW]);
                    if(hammingDist > bigValueL)
                        bigValueL = hammingDist;
                    costL[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += hammingDist;
                }
                else
                    costL[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += bigValueL;

                if(COL+DISP < imgCols)
                {
                    int hammingDist = censusHammingDistance(censusL[imgRows*(COL+DISP) + ROW], censusR[imgRows*COL + ROW]);
                    if(hammingDist > bigValueR)
                        bigValueR = hammingDist;
                    costR[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += hammingDist;
                }
                else
                    costR[imgRows*COL + ROW + imgRows*imgCols*(DISP-1)] += bigValueR;
            }
        }
}

int censusHammingDistance(unsigned char* pxlCensusL, unsigned char* pxlCensusR)
{
    int dist = 0;
    for(int pos = 0; pos < 9*7; ++pos)
        if(pxlCensusL[pos] != pxlCensusR[pos])
            dist++;
    return dist;
}

// Returns those pixels that are not consistent between 2 correspondent stereo disparities
void crossValidation(int* disparity, bool* consistency, int imgRows, int imgCols, int* dispL, int* dispR)
{
    for(int col = 0; col < imgCols; ++col)
        for(int row = 0; row < imgRows; ++row)
        {
            int idxL = imgRows*col + row;
            int idxR = imgRows*(col-dispL[idxL]) + row;
            if(col-dispL[idxL] < 0)
            {
                consistency[idxL] = false;
                disparity[idxL] = -1.0;
            }
            else if(fabs(dispL[idxL] - dispR[idxR]) > imgCols*0.01) // 1% of total width
            {
                consistency[idxL] = false;
                disparity[idxL] = -1.0;
            }
            else
            {
                consistency[idxL] = true;
                disparity[idxL] = dispL[idxL];
            }
        }
}

// Estimates planes from disparity to refine disparity calculations
void computeSegmentPlanes(double* newDisparity, double* planes, int imgRows, int imgCols, int* segmentation, int numSegments, int* disparity, double* scores, bool* consistency, int maxDisparity, int maxIterPlane)
{
    // (a,b,c) per segment = plane.
    for(int idx = 0; idx < numSegments*3; ++idx)
        planes[idx] = 0.0;

    // Counter pixels per segment
    int counter[numSegments];
    // and also gather array of rows, cols and disps per pixel of the segment
    int** idxRow = new int*[numSegments];
    int** idxCol = new int*[numSegments];
    int** dispSeg = new int*[numSegments];
    double** scoresSeg = new double*[numSegments];

    bool hasConsistentPixel[numSegments];

    // Count number of pixels per segment. Discard inconsistencies
    for(int idx = 0; idx < numSegments; ++idx)
        counter[idx] = 0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        if(consistency[idx] == 1.0)
            counter[segmentation[idx]-1]++;

    for(int idx = 0; idx < numSegments; ++idx)
    {
        hasConsistentPixel[idx] = false;
        idxRow[idx] = new int[counter[idx]];
        idxCol[idx] = new int[counter[idx]];
        dispSeg[idx] = new int[counter[idx]];
        scoresSeg[idx] = new double[counter[idx]*maxDisparity];
        counter[idx] = 0;
    }

    // For efficiency: store arrays of ROWS, COLS and DISPARITY VALUES
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = imgRows*COL + ROW;
            if(consistency[idx] == 1.0)
            {
                int SEGMENT = segmentation[idx]-1;

                idxRow[SEGMENT][counter[SEGMENT]] = ROW+1;
                idxCol[SEGMENT][counter[SEGMENT]] = COL+1;
                dispSeg[SEGMENT][counter[SEGMENT]] = disparity[idx];
                for(int DISP = 0; DISP < maxDisparity; ++DISP)
                    scoresSeg[SEGMENT][maxDisparity*counter[SEGMENT] + DISP] = scores[idx + imgRows*imgCols*DISP];
                counter[SEGMENT]++;

                hasConsistentPixel[SEGMENT] = true;
            }
        }

    /*
    mexPrintf("Num. Pixels Segment 302: %d\n", counter[302]);
    for(int i = 0; i < counter[302]; ++i)
    {
        mexPrintf("ROW: %d COL: %d\n", idxRow[302][i], idxCol[302][i]);
        mexPrintf("DISP: %d\n", dispSeg[302][i]);
    }
    */

    // Compute PLANE values per segment
    double* ATxA = new double[9];
    double ATxb[3];
    double x[3];
    for(int SEGMENT = 0; SEGMENT < numSegments; ++SEGMENT)
    {
        double numElems = counter[SEGMENT];

        // If no pixel assigned to current segment, go to next one
        if(numElems == 0)
            continue;

        for(int i = 0; i < 9; ++i)
            ATxA[i] = 0.0;
        for(int i = 0; i < counter[SEGMENT]; ++i)
        {
            ATxA[0] += idxCol[SEGMENT][i]*idxCol[SEGMENT][i];
            ATxA[1] += idxCol[SEGMENT][i]*idxRow[SEGMENT][i];
            ATxA[2] += idxCol[SEGMENT][i];
            ATxA[3] += idxCol[SEGMENT][i]*idxRow[SEGMENT][i];
            ATxA[4] += idxRow[SEGMENT][i]*idxRow[SEGMENT][i];
            ATxA[5] += idxRow[SEGMENT][i];
            ATxA[6] += idxCol[SEGMENT][i];
            ATxA[7] += idxRow[SEGMENT][i];
            ATxA[8]++;
        }
        double detA = det3x3(ATxA);
        inverse3x3(ATxA);

        // Iterate to make a plane free of outliers.
        for(int it = 0; it < maxIterPlane; ++it)
        {
            // Find 'a', 'b' and 'c' by linear least squares.
            // x = (A^T x A)^-1 x (A^T x b)
            if(detA != 0.0)
            {
                ATxb[0] = 0.0; ATxb[1] = 0.0; ATxb[2] = 0.0;
                for(int i = 0; i < counter[SEGMENT]; ++i)
                {
                    ATxb[0] += idxCol[SEGMENT][i]*dispSeg[SEGMENT][i];
                    ATxb[1] += idxRow[SEGMENT][i]*dispSeg[SEGMENT][i];
                    ATxb[2] += dispSeg[SEGMENT][i];
                }
                double* mul = mulMat(3,3,ATxA,3,1,ATxb);
                x[0] = mul[0]; x[1] = mul[1]; x[2] = mul[2];
                delete [] mul;
            }
            else
            {
                x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
                for(int i = 0; i < counter[SEGMENT]; ++i)
                    x[2] += dispSeg[SEGMENT][i];
                x[2] /= counter[SEGMENT];
                break;
            }

            int range = ceil(0.10*maxDisparity); // TODO 0.10 is ok? Too restrictive? Too relaxed?
            // mexPrintf("Max disparity: %d range plane: %d \n", maxDisparity, range);
            // Get the best score within max plane distance.
            for(int ELEM = 0; ELEM < numElems; ++ELEM)
            {
                int planePos = round(x[0]*idxCol[SEGMENT][ELEM] + x[1]*idxRow[SEGMENT][ELEM] + x[2]);
                int limitMin = max(1, planePos - range);
                int limitMax = min(maxDisparity, planePos + range);
                int matchDisp = limitMin;
                double matchScore = scoresSeg[SEGMENT][maxDisparity*ELEM + limitMin-1];
                for(int DISP = limitMin + 1; DISP < limitMax; ++DISP)
                {
                    double newScore = scoresSeg[SEGMENT][maxDisparity*ELEM + DISP - 1];
                    if(newScore < matchScore)
                    {
                        matchScore = newScore;
                        matchDisp = DISP;
                    }
                }
                dispSeg[SEGMENT][ELEM] = matchDisp;
            }
        }

        planes[SEGMENT] = x[0];
        planes[numSegments*1 + SEGMENT] = x[1];
        planes[numSegments*2 + SEGMENT] = x[2];
    }

    // Assign the correspondent values to the other segment pixels (inconsistent)
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = imgRows*COL + ROW;
            int SEGMENT = segmentation[idx]-1;
            if(hasConsistentPixel[SEGMENT])
            {
                //if(consistency[idx] == 1.0)
                //    newDisparity[idx] = disparity[idx];
                //else if(consistency[idx] == -1.0)
                    newDisparity[idx] = (COL+1)*planes[SEGMENT] + (ROW+1)*planes[numSegments*1 + SEGMENT] + planes[numSegments*2 + SEGMENT];

                if(newDisparity[idx] < 1.0)
                    newDisparity[idx] = 1.0;
                else if(newDisparity[idx] > maxDisparity)
                    newDisparity[idx] = maxDisparity;
            }
            else
                newDisparity[idx] = -1.0;
        }

    delete [] ATxA;
    for(int idx = 0; idx < numSegments; ++idx)
    {
        delete [] idxRow[idx];
        delete [] idxCol[idx];
        delete [] dispSeg[idx];
        delete [] scoresSeg[idx];
    }
    delete [] idxRow; delete [] idxCol; delete [] dispSeg; delete [] scoresSeg;
}

// Updates the current disparity by filling inconsistent pixels with closest best matches (first segment, then overall)
void consistentFill(double* newDisparity, int imgRows, int imgCols, double* disparity)
{
    // Pre-processing: Values out of the scope become inconsistent and filled with nearest neighbours
    // - Estimate mean
    double meanDisp = 0.0;
    int count = 0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        if(disparity[idx] >= 1.0)
        {
            meanDisp += disparity[idx];
            count++;
        }
    meanDisp /= count;
    // - Estimate std. deviation
    double devDisp = 0.0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        if(disparity[idx] >= 1.0)
            devDisp += pow(disparity[idx] - meanDisp, 2);
    devDisp = sqrt(devDisp/count);
    // - Outlier out of 2*devDisp are set inconsistent
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        if(disparity[idx] > meanDisp+2.5*devDisp || disparity[idx] < meanDisp-2.5*devDisp)
            disparity[idx] = -1.0;

    // Proceed with filling of inconsistent regions
    for(int row = 0; row < imgRows; ++row)
        for(int col = 0; col < imgCols; ++col)
        {
            int idx = imgRows*col + row;
            if(disparity[idx] != -1.0)
                newDisparity[idx] = disparity[idx];
            else
            {
                //mexPrintf("row: %d col: %d\n", row, col);
                //mexPrintf("Inconsistent pixels found!\n");
                // Check closest valid pixel within the segment
                bool isFound = false;
                int rad = 1;
                int maxRad = (imgCols*imgRows)/2;
                while(!isFound && rad < maxRad)
                {
                    int limitLeft = max(col-rad, 0);
                    int limitRight = min(col+rad, imgCols-1);
                    int limitTop = max(row-rad, 0);
                    int limitBottom = min(row+rad, imgRows-1);

                    // Top line
                    if(row-rad >= 0 && !isFound)
                        for(int radCol = limitLeft; radCol < limitRight; ++radCol)
                        {
                            int newIdx = imgRows*radCol + row;
                            if(disparity[newIdx] != -1.0)
                            {
                                newDisparity[idx] = disparity[newIdx];
                                isFound = true;
                            }
                        }

                    // Bottom line
                    if(row+rad < imgRows && !isFound)
                        for(int radCol = limitLeft; radCol < limitRight; ++radCol)
                        {
                            int newIdx = imgRows*radCol + row;
                            if(disparity[newIdx] != -1.0)
                            {
                                newDisparity[idx] = disparity[newIdx];
                                isFound = true;
                            }
                        }

                    // Left line
                    if(col-rad >= 0 && !isFound)
                        for(int radRow = limitTop; radRow < limitBottom; ++radRow)
                        {
                            int newIdx = imgRows*col + radRow;
                            if(disparity[newIdx] != -1.0)
                            {
                                newDisparity[idx] = disparity[newIdx];
                                isFound = true;
                            }
                        }

                    // Right line
                    if(col+rad < imgCols && !isFound)
                        for(int radRow = limitTop; radRow < limitBottom; ++radRow)
                        {
                            int newIdx = imgRows*col + radRow;
                            if(disparity[newIdx] != -1.0)
                            {
                                newDisparity[idx] = disparity[newIdx];
                                isFound = true;
                            }
                        }
                    rad++;
                }
            }
        }
}

// Estimates similarities among all segments from their planes.
void planeSimilarities(double* segmentSimilarities, int imgRows, int imgCols, int* segmentation, int numSegments, double* planes, double* disparity, bool* consistency, double pImpact)
{
    // PART 1: Estimate planes of previously non consistent segments (which now have been filled)
    int sizeSegment[numSegments];
    for(int idx = 0; idx < numSegments; ++idx)
        sizeSegment[idx] = 0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        sizeSegment[segmentation[idx]-1]++;

    bool noConsistentSegment[numSegments];

    int** idxRow = new int*[numSegments];
    int** idxCol = new int*[numSegments];
    int** dispSeg = new int*[numSegments];
    int counter[numSegments];
    for(int idx = 0; idx < numSegments; ++idx)
    {
        idxRow[idx] = new int[sizeSegment[idx]];
        idxCol[idx] = new int[sizeSegment[idx]];
        dispSeg[idx] = new int[sizeSegment[idx]];
        noConsistentSegment[idx] = false;
        counter[idx] = 0;
    }

    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = imgRows*COL + ROW;
            if(consistency[idx] == 0.0)
            {
                int SEGMENT = (int)segmentation[idx]-1;
                idxRow[SEGMENT][counter[SEGMENT]] = ROW+1;
                idxCol[SEGMENT][counter[SEGMENT]] = COL+1;
                dispSeg[SEGMENT][counter[SEGMENT]] = disparity[idx];
                counter[SEGMENT]++;
                noConsistentSegment[SEGMENT] = true;
            }
        }

    double* ATxA = new double[9];
    double ATxb[3];
    for(int SEGMENT = 0; SEGMENT < numSegments; ++SEGMENT)
        if(noConsistentSegment[SEGMENT] && counter[SEGMENT] != 0)
        {
            for(int i = 0; i < 9; ++i)
                ATxA[i] = 0.0;
            for(int i = 0; i < counter[SEGMENT]; ++i)
            {
                ATxA[0] += idxCol[SEGMENT][i]*idxCol[SEGMENT][i];
                ATxA[1] += idxCol[SEGMENT][i]*idxRow[SEGMENT][i];
                ATxA[2] += idxCol[SEGMENT][i];
                ATxA[3] += idxCol[SEGMENT][i]*idxRow[SEGMENT][i];
                ATxA[4] += idxRow[SEGMENT][i]*idxRow[SEGMENT][i];
                ATxA[5] += idxRow[SEGMENT][i];
                ATxA[6] += idxCol[SEGMENT][i];
                ATxA[7] += idxRow[SEGMENT][i];
                ATxA[8]++;
            }
            double detA = det3x3(ATxA);
            inverse3x3(ATxA);

            // Find 'a', 'b' and 'c' by linear least squares.
            // x = (A^T x A)^-1 x (A^T x b)
            if(detA != 0.0)
            {
                ATxb[0] = 0.0; ATxb[1] = 0.0; ATxb[2] = 0.0;
                for(int i = 0; i < counter[SEGMENT]; ++i)
                {
                    ATxb[0] += idxCol[SEGMENT][i]*dispSeg[SEGMENT][i];
                    ATxb[1] += idxRow[SEGMENT][i]*dispSeg[SEGMENT][i];
                    ATxb[2] += dispSeg[SEGMENT][i];
                }
                double* mul = mulMat(3,3,ATxA,3,1,ATxb);
                planes[SEGMENT] = mul[0]; planes[numSegments*1 + SEGMENT] = mul[1]; planes[numSegments*2 + SEGMENT] = mul[2];
                delete [] mul;
            }
            else
            {
                planes[SEGMENT] = 0.0; planes[numSegments + SEGMENT] = 0.0; planes[numSegments*2 + SEGMENT] = 0.0;
                for(int i = 0; i < counter[SEGMENT]; ++i)
                    planes[numSegments*2 + SEGMENT] += dispSeg[SEGMENT][i];
                planes[numSegments*2 + SEGMENT] /= counter[SEGMENT];
                // mexPrintf("SEGMENT %d; a = %f, b = %f, c = %f\n", SEGMENT, planes[SEGMENT], planes[numSegments + SEGMENT], planes[numSegments*2 + SEGMENT]);
            }
        }

    delete [] ATxA;
    for(int idx = 0; idx < numSegments; ++idx)
    {
        delete [] idxRow[idx];
        delete [] idxCol[idx];
        delete [] dispSeg[idx];
    }
    delete [] idxRow; delete [] idxCol; delete [] dispSeg;

    // Part 2: Computation of average depth for segment
    double meanSegments[numSegments];
    for(int SEGMENT = 0; SEGMENT < numSegments; ++SEGMENT)
        meanSegments[SEGMENT] = 0.0;

    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = imgRows*COL + ROW;
            meanSegments[segmentation[idx]-1] += disparity[idx];
        }

    for(int SEGMENT = 0; SEGMENT < numSegments; ++SEGMENT)
        meanSegments[SEGMENT] /= sizeSegment[SEGMENT];

    // Part 3: Computation of segment similarities/penalties
    double meanDepth = 0.0;
    double meanAngle = 0.0;
    double* pSegment = new double[numSegments*numSegments*2];
    for(int SEG_CURRENT = 0; SEG_CURRENT < numSegments; ++SEG_CURRENT)
    {
        // Get plane parameters.
        double a1 = planes[SEG_CURRENT];
        double b1 = planes[numSegments*1 + SEG_CURRENT];
        double c1 = planes[numSegments*2 + SEG_CURRENT];

        // Access to only once in the comparison (upper triangular matrix) -> avoids redundance
        for(int SEG_COMPARE = SEG_CURRENT+1; SEG_COMPARE < numSegments; ++SEG_COMPARE)
        {
            // Get "plane to compare" parameters.
            double a2 = planes[SEG_COMPARE];
            double b2 = planes[numSegments*1 + SEG_COMPARE];
            double c2 = planes[numSegments*2 + SEG_COMPARE];

            // Get depth-distance penalty.
            double depth = fabs(meanSegments[SEG_CURRENT]- meanSegments[SEG_COMPARE]);
            pSegment[numSegments*SEG_CURRENT + SEG_COMPARE] = depth;
            meanDepth += depth;

            // Get angle penalty (radians).
            double dot = a1*a2 + b1*b2 + c1*c2;
            double cross = sqrt(a1*a1 + b1*b1 + c1*c1)*sqrt(a2*a2 + b2*b2 + c2*c2);
            double angle = acos(dot/cross);
            pSegment[numSegments*SEG_CURRENT + SEG_COMPARE + numSegments*numSegments] = angle;
            meanAngle += angle;
        }
        // mexPrintf("current meanAngle seg %d: %f\n", SEG_CURRENT, meanAngle);
    }

    // Reduce influence of outliers
    int numComparisons = (numSegments*numSegments - numSegments)/2.0;
    meanDepth /= numComparisons;
    meanAngle /= numComparisons;

    double devDepth = 0.0;
    double devAngle = 0.0;
    for(int SEG1 = 0; SEG1 < numSegments; ++SEG1)
        for(int SEG2 = SEG1+1; SEG2 < numSegments; ++SEG2)
        {
            devDepth += pow(pSegment[numSegments*SEG1 + SEG2] - meanDepth, 2.0);
            devAngle += pow(pSegment[numSegments*SEG1 + SEG2 + numSegments*numSegments] - meanAngle, 2.0);
        }

    devDepth = sqrt(devDepth/numComparisons);
    devAngle = sqrt(devAngle/numComparisons);

    double limitDepth = meanDepth + devDepth;
    double limitAngle = meanAngle + devAngle;

    for(int SEG1 = 0; SEG1 < numSegments; ++SEG1)
        for(int SEG2 = SEG1+1; SEG2 < numSegments; ++SEG2)
            if(pSegment[numSegments*SEG1 + SEG2] > limitDepth)
            {
                pSegment[numSegments*SEG1 + SEG2] = limitDepth;
                pSegment[numSegments*SEG1 + SEG2 + numSegments*numSegments] = limitAngle;
            }

    // Normalize [0..pImpact]
    for(int SEG1 = 0; SEG1 < numSegments; ++SEG1)
        for(int SEG2 = SEG1; SEG2 < numSegments; ++SEG2)
        {
            double penalty = 1.0;

            if(SEG1 != SEG2)
            {
                double depthSim = pSegment[numSegments*SEG1 + SEG2];
                double angleSim = pSegment[numSegments*SEG1 + SEG2 + numSegments*numSegments];
                if(limitDepth == 0.0)
                    limitDepth = 1.0;
                if(limitAngle == 0.0)
                    limitAngle = 1.0;
                penalty = exp(-(depthSim/limitDepth) * pImpact); // + angleSim/limitAngle)/2.0

                // If too close to 0 -> 0. No contribution -> wall (dirichlet B.C.).
                if(penalty < pImpact/10.0)
                    penalty = 0.0;
            }
            pSegment[numSegments*SEG1 + SEG2] = penalty;
        }

    // Once penalties among segments are computed... assign among pixels.
    // (Up, Down, Left, Right) - pixels for each pixel in the image.
    for(int COL = 0; COL < imgCols; ++COL)
        for(int ROW = 0; ROW < imgRows; ++ROW)
        {
            int idx = imgRows*COL + ROW;

            if(ROW == 0 || COL == 0 || ROW == imgRows-1 || COL == imgCols-1)
            {
                segmentSimilarities[idx + imgRows*imgCols*0] = 1.0;
                segmentSimilarities[idx + imgRows*imgCols*1] = 1.0;
                segmentSimilarities[idx + imgRows*imgCols*2] = 1.0;
                segmentSimilarities[idx + imgRows*imgCols*3] = 1.0;
                continue;
            }

            // Info current pixel
            int segPxl = segmentation[idx]-1;
            int segTop = segmentation[imgRows*COL + ROW-1]-1;
            int segBottom = segmentation[imgRows*COL + ROW+1]-1;
            int segLeft = segmentation[imgRows*(COL-1) + ROW]-1;
            int segRight = segmentation[imgRows*(COL+1) + ROW]-1;

            // TOP
            if(segPxl < segTop)
                segmentSimilarities[idx + imgRows*imgCols*0] = pSegment[numSegments*segPxl + segTop];
            else
                segmentSimilarities[idx + imgRows*imgCols*0] = pSegment[numSegments*segTop + segPxl];

            // DOWN
            if(segPxl < segBottom)
                segmentSimilarities[idx + imgRows*imgCols*1] = pSegment[numSegments*segPxl + segBottom];
            else
                segmentSimilarities[idx + imgRows*imgCols*1] = pSegment[numSegments*segBottom + segPxl];

            // LEFT
            if(segPxl < segLeft)
                segmentSimilarities[idx + imgRows*imgCols*2] = pSegment[numSegments*segPxl + segLeft];
            else
                segmentSimilarities[idx + imgRows*imgCols*2] = pSegment[numSegments*segLeft + segPxl];

            // RIGHT
            if(segPxl < segRight)
                segmentSimilarities[idx + imgRows*imgCols*3] = pSegment[numSegments*segPxl + segRight];
            else
                segmentSimilarities[idx + imgRows*imgCols*3] = pSegment[numSegments*segRight + segPxl];
        }

   delete [] pSegment;
}


