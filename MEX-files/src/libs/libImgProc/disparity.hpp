#ifndef _DISPARITY_
#define _DISPARITY_

#include <string>

#include "structs.hpp"

void computeLocalDisparities(int* dispL, double* scoresL, int* dispR, double* scoresR, int imgRows, int imgCols, int* segL, int numSegL, int* segR, int numSegR,
                             unsigned char* imageL, unsigned char* imageR, int imgChannels, int maxDisparity, double dispLambda, int baseline, char* mCost);

template <class T>
void calculateIndividualCosts(int imgRows, int imgCols, int maxDisparity, int channel, T* imageL, double* costL, T* imageR, double* costR);

void imageCensus(unsigned char** census, int imgRows, int imgCols, unsigned char* image);
void censusCost(int imgRows, int imgCols, int maxDisparity, unsigned char** censusL, double* matchCostL, unsigned char** censusR, double* matchCostR);
int censusHammingDistance(unsigned char* pxlCensusL, unsigned char* pxlCensusR);

void crossValidation(int* disparity, bool* consistency, int imgRows, int imgCols, int* dispL, int* dispR);

void computeSegmentPlanes(double* newDisparity, double* planes, int imgRows, int imgCols, int* segmentation, int numSegments, int* disparity, double* scores, bool* consistency, int maxDisparity, int maxIterPlane);

void consistentFill(double* newDisparity, int imgRows, int imgCols, double* disparity);

void planeSimilarities(double* segmentSimilarities, int imgRows, int imgCols, int* segmentation, int numSegments, double* planes, double* disparity, bool* consistency, double pImpact);

#endif

