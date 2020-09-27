#include <math.h>
#include <vector>
#include <ctime>
#include <set>
#include <queue>
#include <limits>
#include <algorithm>
#include <iostream>

#include "mex.h"

#include "differentiation.hpp"
#include "segmentation.hpp"

using namespace std;

// Returns a segmented image with K-Means algorithm with 7D/9D data.
void KMeans9D(int* segmentation, int imgRows, int imgCols, unsigned char* imageIn, int imgChannels, double* zIn, double* uIn, double* vIn, double* wIn, InputKMeans input)
{
    // Dimensions of dimensions to cluster
    int numDims = imgChannels + 6;

    // Pre-rocessing: normalise input information
    double *image = new double[imgRows*imgCols*imgChannels];
    for(int idx = 0; idx < imgCols*imgRows*imgChannels; ++idx)
        image[idx] = (double)imageIn[idx] / 255.0;

    double* X = new double[imgRows*imgCols];
    double* Y = new double[imgRows*imgCols];
    double* Z = new double[imgRows*imgCols];
    double* u = new double[imgRows*imgCols];
    double* v = new double[imgRows*imgCols];
    double* w = new double[imgRows*imgCols];
    for(int col = 0; col < imgCols; ++col)
        for(int row = 0; row < imgRows; ++row)
        {
            int idx = imgRows*col + row;
            X[idx] = col;
            Y[idx] = row;
        }

    double XMin = *min_element(X,X+imgRows*imgCols); double XMax = *max_element(X,X+imgRows*imgCols);
    double XDelta = (XMax - XMin);
    if(XDelta == 0.0)
        XDelta = 1.0;
    double YMin = *min_element(Y,Y+imgRows*imgCols); double YMax = *max_element(Y,Y+imgRows*imgCols);
    double YDelta = (YMax - YMin);
    if(YDelta == 0.0)
        YDelta = 1.0;
    double ZMin = *min_element(zIn,zIn+imgRows*imgCols); double ZMax = *max_element(zIn,zIn+imgRows*imgCols);
    double ZDelta = (ZMax - ZMin);
    if(ZDelta == 0.0)
        ZDelta = 1.0;
    double uMin = *min_element(uIn,uIn+imgRows*imgCols); double uMax = *max_element(uIn,uIn+imgRows*imgCols);
    double uDelta = (uMax - uMin);
    if(uDelta == 0.0)
        uDelta = 1.0;
    double vMin = *min_element(vIn,vIn+imgRows*imgCols); double vMax = *max_element(vIn,vIn+imgRows*imgCols);
    double vDelta = (vMax - vMin);
    if(vDelta == 0.0)
        vDelta = 1.0;
    double wMin = *min_element(wIn,wIn+imgRows*imgCols); double wMax = *max_element(wIn,wIn+imgRows*imgCols);
    double wDelta = (wMax - wMin);
    if(wDelta == 0.0)
        wDelta = 1.0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
    {
        X[idx] = (X[idx] - XMin)/XDelta;
        Y[idx] = (Y[idx] - YMin)/YDelta;
        Z[idx] = (zIn[idx] - ZMin)/ZDelta;
        u[idx] = (uIn[idx] - uMin)/uDelta;
        v[idx] = (vIn[idx] - vMin)/vDelta;
        w[idx] = (wIn[idx] - wMin)/wDelta;
    }

    double ratio = (double)imgRows/(double)imgCols;
    double centroidsXrow = round(ratio*sqrt(input.numK/ratio));
    double centroidsXcol = round(sqrt(input.numK/ratio));
    if(centroidsXrow * centroidsXcol > input.numK)
        centroidsXrow--;

    double proportionRow = 1.0/centroidsXrow;
    double proportionCol = 1.0/centroidsXcol;
    // mexPrintf("Max. number of initial centroids in K-means: %d\n", input.numK);
    double centroids[input.numK*numDims];

    // Posterior gradient purposes
    double* img_xy = new double[imgRows*imgCols*2];
    gradient(imgRows, imgCols, imageIn, img_xy);

    float segmentRow = proportionRow*imgRows;
    float segmentCol = proportionCol*imgCols;

    int K = 0;
    for(float col = segmentCol/2; col < imgCols; col += segmentCol)
        for(float row = segmentRow/2; row < imgRows; row += segmentRow)
        {
            int intRow = (int)row;
            int intCol = (int)col;

            // Check gradient in 3x3 windows and assign lowest gradient pixel:
            double minGrad = img_xy[imgRows*intCol+ intRow];
            int auxRow = intRow; int auxCol = intCol;
            for(int rowSeg = intRow-1; rowSeg <= intRow+1; ++rowSeg)
                for(int colSeg = intCol-1; colSeg <= intCol+1; ++colSeg)
                {
                    if(img_xy[imgRows*rowSeg + colSeg] < minGrad)
                    {
                        auxRow = rowSeg;
                        auxCol = colSeg;
                        minGrad = img_xy[imgRows*rowSeg + colSeg];
                    }
                }
            intRow = auxRow; intCol = auxCol;

            int idx = imgRows*intCol+ intRow;
            // Structure of dimensional space to cluster --> 7D: [X Y Z u v w i] or 9D: [X Y Z u v w r g b]
            centroids[numDims*K + 0] = X[idx]; centroids[numDims*K + 1] = Y[idx]; centroids[numDims*K + 2] = Z[idx];
            centroids[numDims*K + 3] = u[idx]; centroids[numDims*K + 4] = v[idx]; centroids[numDims*K + 5] = w[idx];
            centroids[numDims*K + 6] = image[idx];
            if(numDims > 7)
            {
                centroids[numDims*K + 7] = image[idx + imgRows*imgCols];
                centroids[numDims*K + 8] = image[idx + imgRows*imgCols*2];
            }

            K++;
        }

    // Actual number of centroids (rounded -> sqrt (int))
    input.numK = K;
    // mexPrintf("Actual number of centroids: %d\n", input.numK);

    // Counter for number of pixels assigned to a centroid
    int counter[input.numK];
    // Run K-means.
    for (int it = 0; it < input.numIterKMeans; ++it)
    {
        // Assign a centroid per pixel.
        for(int COL = 0; COL < imgCols; ++COL)
            for(int ROW = 0; ROW < imgRows; ++ROW)
            {
                int idx = imgRows*COL + ROW;

                double bestSum = numeric_limits<double>::infinity();
                int bestCentroid = -1;
                for(int CENTROID = 0; CENTROID < input.numK; ++CENTROID)
                {
                    // Discard centroids out of 2*windowsize distance in the image.
                    if(fabs(X[idx] - centroids[numDims*CENTROID + 0]) > 3*proportionCol && fabs(Y[idx] - centroids[numDims*CENTROID + 1]) > 3*proportionRow)
                        continue;

                    double sum = 0.0;
                    sum += fabs(X[idx] - centroids[numDims*CENTROID + 0]) * 0.25;
                    sum += fabs(Y[idx] - centroids[numDims*CENTROID + 1]) * 0.25;
                    sum += fabs(Z[idx] - centroids[numDims*CENTROID + 2]);
                    sum += fabs(u[idx] - centroids[numDims*CENTROID + 3]);
                    sum += fabs(v[idx] - centroids[numDims*CENTROID + 4]);
                    sum += fabs(w[idx] - centroids[numDims*CENTROID + 5]);
                    // Treatment whether 7 or 9 dimension feature vectors.
                    if(numDims == 7)
                    {
                        sum += fabs(image[idx] - centroids[numDims*CENTROID + 6]);
                    }
                    else if(numDims > 7)
                    {
                        sum += 0.33*fabs(image[idx] - centroids[numDims*CENTROID + 6]);
                        sum += 0.33*fabs(image[idx + imgRows*imgCols] - centroids[numDims*CENTROID + 7]);
                        sum += 0.33*fabs(image[idx + imgRows*imgCols*2] - centroids[numDims*CENTROID + 8]);
                    }

                    if(sum < bestSum)
                    {
                        bestSum = sum;
                        bestCentroid = CENTROID;
                    }
                }
                // mexPrintf("%d\n", bestCentroid);
                segmentation[idx] = bestCentroid;
            }

        for(int count = 0; count < input.numK*numDims; ++count)
            centroids[count] = 0;
        for(int count = 0; count < input.numK; ++count)
            counter[count] = 0;
        // Recalculate centroids
        for(int COL = 0; COL < imgRows; ++COL)
            for(int ROW = 0; ROW < imgRows; ++ROW)
            {
                int idx = imgRows*COL + ROW;
                centroids[numDims*segmentation[idx] + 0] += X[idx];
                centroids[numDims*segmentation[idx] + 1] += Y[idx];
                centroids[numDims*segmentation[idx] + 2] += Z[idx];
                centroids[numDims*segmentation[idx] + 3] += u[idx];
                centroids[numDims*segmentation[idx] + 4] += v[idx];
                centroids[numDims*segmentation[idx] + 5] += w[idx];
                centroids[numDims*segmentation[idx] + 6] += image[idx];
                if(numDims > 7)
                {
                    centroids[numDims*segmentation[idx] + 7] += image[idx + imgRows*imgCols];
                    centroids[numDims*segmentation[idx] + 8] += image[idx + imgRows*imgCols*2];
                }
                // Add to the number of pixels belonging to the centroid.
                counter[segmentation[idx]]++;
            }

        // Average values
        for(int CENTROID = 0; CENTROID < input.numK; ++CENTROID)
            for(int DIM = 0; DIM < numDims; ++DIM)
                centroids[numDims*CENTROID + DIM] /= counter[CENTROID];
    }

    // Minimum segment = 1 and not 0.
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        segmentation[idx]++;

    delete [] img_xy;
    delete [] image;
    delete [] X; delete [] Y; delete [] Z;
    delete [] u; delete [] v; delete [] w;
}

// Non-grouped segmentation transformed in a connected component segmented image
void superPixelSegmentation(int imgRows, int imgCols, int imgChannels, int* seg, unsigned char* imgIn, double* zIn, double* uIn, double* vIn, double* wIn)
{
    // Pre-rocessing: normalise input information
    double *img = new double[imgRows*imgCols*imgChannels];
    for(int idx = 0; idx < imgCols*imgRows*imgChannels; ++idx)
        img[idx] = (double)imgIn[idx] / 255.0;
    double* Z = new double[imgRows*imgCols];
    double* u = new double[imgRows*imgCols];
    double* v = new double[imgRows*imgCols];
    double* w = new double[imgRows*imgCols];
    double ZMin = *min_element(zIn,zIn+imgRows*imgCols); double ZMax = *max_element(zIn,zIn+imgRows*imgCols);
    double ZDelta = (ZMax - ZMin);
    if(ZDelta == 0.0)
        ZDelta = 1.0;
    double uMin = *min_element(uIn,uIn+imgRows*imgCols); double uMax = *max_element(uIn,uIn+imgRows*imgCols);
    double uDelta = (uMax - uMin);
    if(uDelta == 0.0)
        uDelta = 1.0;
    double vMin = *min_element(vIn,vIn+imgRows*imgCols); double vMax = *max_element(vIn,vIn+imgRows*imgCols);
    double vDelta = (vMax - vMin);
    if(vDelta == 0.0)
        vDelta = 1.0;
    double wMin = *min_element(wIn,wIn+imgRows*imgCols); double wMax = *max_element(wIn,wIn+imgRows*imgCols);
    double wDelta = (wMax - wMin);
    if(wDelta == 0.0)
        wDelta = 1.0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
    {
        Z[idx] = (zIn[idx] - ZMin)/ZDelta;
        u[idx] = (uIn[idx] - uMin)/uDelta;
        v[idx] = (vIn[idx] - vMin)/vDelta;
        w[idx] = (wIn[idx] - wMin)/wDelta;
    }

    // Once normalised to [0..1] discretise in N number of bins -> test: nBins = 8
    // int nBins = 32;
    /*
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
    {
        Z[idx] = floor(Z[idx] * nBins*0.99);
        u[idx] = floor(u[idx] * nBins*0.99);
        v[idx] = floor(v[idx] * nBins*0.99);
        w[idx] = floor(w[idx] * nBins*0.99);
    }
    for(int idx = 0; idx < imgCols*imgRows*imgChannels; ++idx)
        img[idx] = floor(img[idx] * nBins*0.99);
    */

    int imgLabels = twoPassConnectedComponentAlgorithm(imgRows, imgCols, seg);

    //Threshold for considering a segment too small
    int thSize = max((int)pow(min(imgRows,imgCols)*0.01, 2.0), 9);
    //mexPrintf("Connected labels: %d\n", imgLabels);
    //mexPrintf("Minimum number of pixels per segment: %d\n", thSize);

    agglomerativeClustering(imgRows, imgCols, imgChannels, imgLabels, seg, img, Z, u, v, w, thSize, false);

    delete [] img; delete [] Z; delete [] u; delete [] v; delete [] w;
}

int twoPassConnectedComponentAlgorithm(int imgRows, int imgCols, int* seg)
{
    vector<int> labels(imgRows*imgCols);
    for(int i = 0; i < imgRows*imgCols; ++i)
        labels[i] = -1;
    vector< set<int> > connectedNeighbours;
    vector< set<int> > nonConnectedNeighbours;

    // First pass
    int iLabel = 0;
    for(int row = 0; row < imgRows; ++row)
        for(int col = 0; col < imgCols; ++col)
        {
            int idx = imgRows*col + row;
            int idxLeft = imgRows*(col-1) + row;
            int idxTop = imgRows*col + row-1;

            bool sameLeft = (col-1 >= 0);
            if(sameLeft)
                sameLeft = (seg[idx] == seg[idxLeft]);
            bool sameTop = (row-1 >= 0);
            if(sameTop)
                sameTop = (seg[idx] == seg[idxTop]);

            // Check left pxl if same segment/label
            if(sameLeft)
                labels[idx] = labels[idxLeft];

            // Check top pxl if same segment/label
            if(sameTop)
            {
                bool assignedTop = false;
                if(!sameLeft || labels[idxTop] < labels[idxLeft])
                {
                    labels[idx] = labels[idxTop];
                    assignedTop = true;
                }

                // Store connectivy in case both neighbours have the same label
                if(assignedTop && sameLeft)
                {
                    connectedNeighbours[labels[idx]].insert(labels[idxLeft]);
                    connectedNeighbours[labels[idxLeft]].insert(labels[idx]);
                }
                else if(sameLeft)
                {
                    connectedNeighbours[labels[idx]].insert(labels[idxTop]);
                    connectedNeighbours[labels[idxTop]].insert(labels[idx]);
                }
            }

            // Assign new label
            if(!sameLeft && !sameTop)
            {
                // New label
                labels[idx] = iLabel;

                // Add new label in connected/nonConnected neighbours
                set<int> newSetConnected;
                set<int> newSetNonConnected;
                connectedNeighbours.push_back(newSetConnected);
                connectedNeighbours[iLabel].insert(iLabel);
                nonConnectedNeighbours.push_back(newSetNonConnected);

                iLabel++;
            }
        }

    // Second pass
    vector<int> counterLabel(iLabel);
    vector<int> updateLabel(iLabel);
    for(int i = 0; i < iLabel; ++i)
    {
        counterLabel[i] = 0;
        updateLabel[i] = -1;
    }
    for(int row = 0; row < imgRows; ++row)
        for(int col = 0; col < imgCols; ++col)
        {
            int idx = imgRows*col + row;
            int assignedLabel;

            if(updateLabel[labels[idx]] == -1)
            {
                // Go through all labels and first one which is connected (lowest one) assign to final segmentation.
                vector<int> search(connectedNeighbours[labels[idx]].begin(), connectedNeighbours[labels[idx]].end());
                vector<int>::iterator it;
                set<int> candidates;

                while(!search.empty())
                {
                    // mexPrintf("size search: %d, size candidates: %d \n", search.size(), candidates.size());
                    int currentElem = *search.begin();
                    // mexPrintf("currentElem: %d \n", currentElem);
                    search.erase(search.begin());
                    if(candidates.find(currentElem) == candidates.end())
                    {
                        // Add to the candidates list
                        candidates.insert(currentElem);
                        // mexPrintf("new candidate: %d\n", currentElem);
                        it = search.end();
                        search.insert(it, connectedNeighbours[currentElem].begin(), connectedNeighbours[currentElem].end());
                    }
                }
                assignedLabel = *candidates.begin();
                // Update the label it assigns
                updateLabel[labels[idx]] = assignedLabel;
                // make the connected labels of the current one only the assigned one (should be slightly more efficient this step)
                connectedNeighbours[labels[idx]].empty();
                connectedNeighbours[labels[idx]].insert(assignedLabel);
            }
            else
            {
                assignedLabel = updateLabel[labels[idx]];
            }

            labels[idx] = assignedLabel;
            counterLabel[assignedLabel]++;
        }

    // Remove empty segments and label with lowest numbers
    for(int i = 0; i < iLabel; ++i)
        updateLabel[i] = -1;

    // Reassign from [0..maxRealLabel]
    int newLabel = 0;
    for(int i = 0; i < iLabel; ++i)
    {
        if(counterLabel[i] > 0)
        {
            updateLabel[i] = newLabel;
            newLabel++;
        }
    }
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        labels[idx] = updateLabel[labels[idx]];

    // Assign labels to final connected segmentation
    for(int i = 0; i < imgRows*imgCols; ++i)
        seg[i] = labels[i];

    return newLabel;
}

void agglomerativeClustering(int imgRows, int imgCols, int imgChannels, int imgLabels, int* seg,
                             double* img, double* Z, double* u, double* v, double* w, int thSize, bool similarityAwareness, double mergeTh)
{
    // Auxiliar labeling result
    vector<int> labels(imgRows*imgCols);
    for(int i = 0; i < imgRows*imgCols; ++i)
        labels[i] = seg[i];

    // Declare features of each segment
    // NOTE: right now only intensity/colour differences
    // NOTE #2: preMean, mean and var are used for merging purposes on (similarityAwareness on)
    int numFeatures = imgChannels + 4;
    vector< vector<double> > features;
    features.resize(2*imgLabels);
    vector< vector<double> > mean;
    mean.resize(2*imgLabels);
    vector< vector<double> > var;
    var.resize(2*imgLabels);
    for (int i = 0; i < 2*imgLabels; ++i)
    {
        features[i].resize(numFeatures);
        mean[i].resize(numFeatures);
        var[i].resize(numFeatures);
    }

    // Nonconnected components
    vector< set<int> > nonConnectedNeighbours(imgLabels);
    for(int i = 0; i < nonConnectedNeighbours.size(); ++i)
        nonConnectedNeighbours[i].clear();
    for(int col = 0; col < imgCols; ++col)
        for(int row = 0; row < imgRows; ++row)
        {
            int idx = imgRows*col + row;
            int idxLeft = imgRows*(col-1) + row;
            int idxTop = imgRows*col + row-1;

            if(col-1 >= 0)
                if(labels[idx] != labels[idxLeft])
                {
                    nonConnectedNeighbours[labels[idx]].insert(labels[idxLeft]);
                    nonConnectedNeighbours[labels[idxLeft]].insert(labels[idx]);
                }

            if(row-1 >= 0)
                if(labels[idx] != labels[idxTop])
                {
                    nonConnectedNeighbours[labels[idx]].insert(labels[idxTop]);
                    nonConnectedNeighbours[labels[idxTop]].insert(labels[idx]);
                }
        }

    // For efficiency, store the pixels that are assigned to a given segment
    vector< set<int> > listPxlSegment(imgLabels);
    for(int label = 0; label < imgLabels; ++label)
        listPxlSegment[label].clear();
    for(int idx = 0; idx < imgCols*imgRows; ++idx)
        listPxlSegment[labels[idx]].insert(idx);

    // Define features per segment
    for(int i = 0; i < 2*imgLabels; ++i)
        for(int f = 0; f < numFeatures; ++f)
        {
            features[i][f] = 0.0;
            mean[i][f] = 0.0;
            var[i][f] = 0.0;
        }

    for(int row = 0; row < imgRows; ++row)
        for(int col = 0; col < imgCols; ++col)
        {
            int idx = imgRows*col + row;
            features[labels[idx]][0] += Z[idx];
            features[labels[idx]][1] += u[idx];
            features[labels[idx]][2] += v[idx];
            features[labels[idx]][3] += w[idx];
            for(int c = 0; c < imgChannels; ++c)
                features[labels[idx]][4+c] += img[idx];
        }

    for(int i = 0; i < imgLabels; ++i)
    {
        int size = listPxlSegment[i].size();
        if(size > 0)
        {
            features[i][0] /= size;
            features[i][1] /= size;
            features[i][2] /= size;
            features[i][3] /= size;
            for(int c = 0; c < imgChannels; ++c)
                features[i][4+c] /= size;
        }
    }

    // Only use std deviation to ignore outliers in the mean if (similarityAwareness)
    if(similarityAwareness)
    {
        // Variance estimation
        double diff;
        for(int i = 0; i < imgRows*imgCols; ++i)
        {
            diff = Z[i] - features[labels[i]][0];
            var[labels[i]][0] += diff*diff;
            diff = u[i] - features[labels[i]][1];
            var[labels[i]][1] += diff*diff;
            diff = v[i] - features[labels[i]][2];
            var[labels[i]][2] += diff*diff;
            diff = w[i] - features[labels[i]][3];
            var[labels[i]][3] += diff*diff;
            for(int c = 0; c < imgChannels; ++c)
            {
                diff = img[i + imgCols*imgRows*c] - features[labels[i]][4+c];
                var[labels[i]][4+c] += diff*diff;
            }
        }
        for(int i = 0; i < imgLabels; ++i)
        {
            int sizeSeg = listPxlSegment[i].size();
            if(sizeSeg > 0)
            {
                var[i][0] /= sizeSeg; var[i][0] = sqrt(var[i][0]);
                var[i][1] /= sizeSeg; var[i][1] = sqrt(var[i][1]);
                var[i][2] /= sizeSeg; var[i][2] = sqrt(var[i][2]);
                var[i][3] /= sizeSeg; var[i][3] = sqrt(var[i][3]);
                for(int c = 0; c < imgChannels; ++c)
                {
                    var[i][4+c] /= sizeSeg;
                    var[i][4+c] = sqrt(var[i][4+c]);
                }
            }
        }

        // Re-estimate the mean excluding those points out of the variance (one-pass outlier removal)
        vector< vector<int> > size;
        size.resize(imgLabels);
        for (int i = 0; i < imgLabels; ++i)
        {
            size[i].resize(numFeatures);
            for (int f = 0; f < numFeatures; ++f)
                size[i][f] = 0;
        }
        for(int i = 0; i < imgRows*imgCols; ++i)
        {
            if(fabs(features[labels[i]][0] - Z[i]) <= var[labels[i]][0])
            {
                mean[labels[i]][0] += Z[i];
                size[labels[i]][0]++;
            }
            if(fabs(features[labels[i]][1] - u[i]) <= var[labels[i]][1])
            {
                mean[labels[i]][1] += u[i];
                size[labels[i]][1]++;
            }
            if(fabs(features[labels[i]][2] - v[i]) <= var[labels[i]][2])
            {
                mean[labels[i]][2] += v[i];
                size[labels[i]][2]++;
            }
            if(fabs(features[labels[i]][3] - w[i]) <= var[labels[i]][3])
            {
                mean[labels[i]][3] += w[i];
                size[labels[i]][3]++;
            }
            for(int c = 0; c < imgChannels; ++c)
            {
                if(fabs(features[labels[i]][4+c] - img[i + imgCols*imgRows*c]) <= var[labels[i]][4+c])
                {
                    mean[labels[i]][4+c] += img[i + imgCols*imgRows*c];
                    size[labels[i]][4+c]++;
                }
            }
        }

        for(int i = 0; i < imgLabels; ++i)
            for(int f = 0; f < numFeatures; ++f)
                features[i][f] = 0.0;
        for(int i = 0; i < imgLabels; ++i)
        {
            features[i][0] = mean[i][0] / size[i][0];
            features[i][1] = mean[i][1] / size[i][1];
            features[i][2] = mean[i][2] / size[i][2];
            features[i][3] = mean[i][3] / size[i][3];
            for(int c = 0; c < imgChannels; ++c)
                features[i][4+c] = mean[i][4+c] / size[i][4+c];
        }
    }

    // Priority queue of distance comparison (less distance more priority)
    priority_queue<Neighbours, vector<Neighbours>, distComparison> priorityComparison;
    for(int i = 0; i < imgLabels; ++i)
        if((listPxlSegment[i].size() < thSize && listPxlSegment[i].size() > 0) || similarityAwareness)
        {
            // Check all non connected neighbouring labels and merge with biggest one (with most similar feature vector)
            vector<int> search(nonConnectedNeighbours[i].begin(), nonConnectedNeighbours[i].end());
            for(vector<int>::iterator it = search.begin(); it != search.end(); ++it)
            {
                int neighbour = *it;

                // Compute similarity distance between neighbours
                double auxDistanceFeature = 0.0;
                for(int f = 0; f < numFeatures; ++f)
                    auxDistanceFeature += fabs(features[i][f] - features[neighbour][f]);

                Neighbours n;
                n.n1 = i;
                n.n2 = neighbour;
                n.dist = auxDistanceFeature;
                priorityComparison.push(n);
            }
    }

    // Comparisons between all neighbours inserted into a priority queue
    vector<bool> validComparison(2*imgLabels);
    for(int i = 0; i < imgLabels; ++i)
        validComparison[i] = true;
    for(int i = imgLabels; i < 2*imgLabels; ++i)
        validComparison[i] = false;

    // Iterative until priority queue with comparisons is empty, as it would mean that no more mergings need to be done.
    int numIter = 0;
    while(!priorityComparison.empty())
    {
        // Merge clusters with highest priority
        Neighbours n = priorityComparison.top();
        priorityComparison.pop();

        // Ignore non-valid comparisons (already merged)
        if(!validComparison[n.n1] || !validComparison[n.n2])
            continue;

        // If we merge by similarity, ignore distances above the given threshold
        if(n.dist >= mergeTh)
            continue;
        //else if(similarityAwareness)
        //    cout << "MERGES! " << n.n1 << " with " << n.n2 << endl;

        validComparison[n.n1] = false;
        validComparison[n.n2] = false;
        validComparison[imgLabels] = true;

        // Update features of the new segment
        int numN1 = listPxlSegment[n.n1].size();
        int numN2 = listPxlSegment[n.n2].size();
        features[imgLabels][0] = (features[n.n1][0] * numN1 + features[n.n2][0] * numN2) / (numN1+numN2);
        features[imgLabels][1] = (features[n.n1][1] * numN1 + features[n.n2][1] * numN2) / (numN1+numN2);
        features[imgLabels][2] = (features[n.n1][2] * numN1 + features[n.n2][2] * numN2) / (numN1+numN2);
        features[imgLabels][3] = (features[n.n1][3] * numN1 + features[n.n2][3] * numN2) / (numN1+numN2);
        for(int c = 0; c < imgChannels; ++c)
            features[imgLabels][4+c] = (features[n.n1][4+c] * numN1 + features[n.n2][4+c] * numN2) / (numN1+numN2);

        // Update set of pxls per segment
        set<int> newSetPxls;
        listPxlSegment.push_back(newSetPxls);
        listPxlSegment[imgLabels].insert(listPxlSegment[n.n1].begin(), listPxlSegment[n.n1].end());
        listPxlSegment[imgLabels].insert(listPxlSegment[n.n2].begin(), listPxlSegment[n.n2].end());
        listPxlSegment[n.n1].clear();
        listPxlSegment[n.n2].clear();

        // Update set of neighbours per segment
        set<int> newSetNeighbour;
        nonConnectedNeighbours.push_back(newSetNeighbour);
        for(set<int>::iterator itNeighbour = nonConnectedNeighbours[n.n1].begin(); itNeighbour != nonConnectedNeighbours[n.n1].end(); ++itNeighbour)
        {
            if(validComparison[*itNeighbour])
            {
                nonConnectedNeighbours[imgLabels].insert(*itNeighbour);
                nonConnectedNeighbours[*itNeighbour].insert(imgLabels);
            }
        }
        for(set<int>::iterator itNeighbour = nonConnectedNeighbours[n.n2].begin(); itNeighbour != nonConnectedNeighbours[n.n2].end(); ++itNeighbour)
        {
            if(validComparison[*itNeighbour])
            {
                nonConnectedNeighbours[imgLabels].insert(*itNeighbour);
                nonConnectedNeighbours[*itNeighbour].insert(imgLabels);
            }
        }
        nonConnectedNeighbours[n.n1].clear();
        nonConnectedNeighbours[n.n2].clear();

        // Only use std deviation to ignore outliers in the mean if (similarityAwareness)
        if(similarityAwareness)
        {
            // Variance estimation
            for(set<int>::iterator itNeighbour = listPxlSegment[imgLabels].begin(); itNeighbour != listPxlSegment[imgLabels].end(); ++itNeighbour)
            {
                int pxl = *itNeighbour;
                double diff;
                diff = Z[pxl] - features[imgLabels][0];
                var[imgLabels][0] += diff*diff;
                diff = u[pxl] - features[imgLabels][1];
                var[imgLabels][1] += diff*diff;
                diff = v[pxl] - features[imgLabels][2];
                var[imgLabels][2] += diff*diff;
                diff = w[pxl] - features[imgLabels][3];
                var[imgLabels][3] += diff*diff;
                for(int c = 0; c < imgChannels; ++c)
                {
                    diff = img[pxl + imgCols*imgRows*c] - features[imgLabels][4+c];
                    var[imgLabels][4+c] += diff*diff;
                }
            }

            int sizeSeg = listPxlSegment[imgLabels].size();
            {
                var[imgLabels][0] /= sizeSeg; var[imgLabels][0] = sqrt(var[imgLabels][0]);
                var[imgLabels][1] /= sizeSeg; var[imgLabels][1] = sqrt(var[imgLabels][1]);
                var[imgLabels][2] /= sizeSeg; var[imgLabels][2] = sqrt(var[imgLabels][2]);
                var[imgLabels][3] /= sizeSeg; var[imgLabels][3] = sqrt(var[imgLabels][3]);
                for(int c = 0; c < imgChannels; ++c)
                {
                    var[imgLabels][4+c] /= sizeSeg;
                    var[imgLabels][4+c] = sqrt(var[imgLabels][4+c]);
                }
            }

            // Re-estimate the mean excluding those points out of the variance (one-pass outlier removal)
            vector<int> size;
            size.resize(numFeatures);
            for (int f = 0; f < numFeatures; ++f)
                size[f] = 0;
            for(set<int>::iterator itNeighbour = listPxlSegment[imgLabels].begin(); itNeighbour != listPxlSegment[imgLabels].end(); ++itNeighbour)
            {
                int pxl = *itNeighbour;
                if(fabs(features[imgLabels][0] - Z[pxl]) <= var[imgLabels][0])
                {
                    mean[imgLabels][0] += Z[pxl];
                    size[0]++;
                }
                if(fabs(features[imgLabels][1] - u[pxl]) <= var[imgLabels][1])
                {
                    mean[imgLabels][1] += u[pxl];
                    size[1]++;
                }
                if(fabs(features[imgLabels][2] - v[pxl]) <= var[imgLabels][2])
                {
                    mean[imgLabels][2] += v[pxl];
                    size[2]++;
                }
                if(fabs(features[imgLabels][3] - w[pxl]) <= var[imgLabels][3])
                {
                    mean[imgLabels][3] += w[pxl];
                    size[3]++;
                }
                for(int c = 0; c < imgChannels; ++c)
                {
                    if(fabs(features[imgLabels][4+c] - img[pxl + imgCols*imgRows*c]) <= var[imgLabels][4+c])
                    {
                        mean[imgLabels][4+c] += img[pxl + imgCols*imgRows*c];
                        size[4+c]++;
                    }
                }
            }

            for(int f = 0; f < numFeatures; ++f)
                features[imgLabels][f] = 0.0;
            features[imgLabels][0] = mean[imgLabels][0] / size[0];
            features[imgLabels][1] = mean[imgLabels][1] / size[1];
            features[imgLabels][2] = mean[imgLabels][2] / size[2];
            features[imgLabels][3] = mean[imgLabels][3] / size[3];
            for(int c = 0; c < imgChannels; ++c)
                features[imgLabels][4+c] = mean[imgLabels][4+c] / size[4+c];
        }

        // Update distances
        for(set<int>::iterator itNeighbour = nonConnectedNeighbours[imgLabels].begin(); itNeighbour != nonConnectedNeighbours[imgLabels].end(); ++itNeighbour)
        {
            int neighbour = *itNeighbour;

            // Ignore comparison with a segment which has already been merged or both segments are big enough
            if(!validComparison[neighbour] || (listPxlSegment[neighbour].size() >= thSize && listPxlSegment[imgLabels].size() >= thSize &&  !similarityAwareness))
                continue;

            // Compute similarity distance between neighbours
            double auxDistanceFeature = 0.0;
            for(int f = 0; f < numFeatures; ++f)
                auxDistanceFeature += fabs(features[imgLabels][f] - features[neighbour][f]);

            Neighbours nAdd;
            nAdd.n1 = neighbour;
            nAdd.n2 = imgLabels;
            nAdd.dist = auxDistanceFeature;
            priorityComparison.push(nAdd);
        }

        imgLabels++;
        numIter++;
    }

    // Assign pixels to its correspondent segment
    for(int segIdx = 0; segIdx < listPxlSegment.size(); ++segIdx)
        for(set<int>::iterator pxlIdx = listPxlSegment[segIdx].begin(); pxlIdx != listPxlSegment[segIdx].end(); ++pxlIdx)
            labels[*pxlIdx] = segIdx;

    // Remove empty segments and label with lowest numbers
    vector<int> updateLabel(imgLabels);
    for(int i = 0; i < imgLabels; ++i)
        updateLabel[i] = -1;

    // Reassign from [0..maxRealLabel]
    int newLabel = 0;
    for(int i = 0; i < imgLabels; ++i)
    {
        if(listPxlSegment[i].size() > 0)
        {
            updateLabel[i] = newLabel;
            newLabel++;
        }
    }
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
        labels[idx] = updateLabel[labels[idx]];

    // Assign labels to final connected segmentation
    for(int i = 0; i < imgRows*imgCols; ++i)
        seg[i] = labels[i];
}

// Merges all adjacent segments that are very similar giving a minimum similarity threshold value
void mergeSegments(int imgRows, int imgCols, int imgChannels, int imgLabels, int* segmentation, unsigned char* imgIn, double* zIn, double* uIn, double* vIn, double* wIn, double mergeTh)
{
    // Pre-rocessing: normalise input information
    double *img = new double[imgRows*imgCols*imgChannels];
    for(int idx = 0; idx < imgCols*imgRows*imgChannels; ++idx)
        img[idx] = (double)imgIn[idx] / 255.0;
    double* Z = new double[imgRows*imgCols];
    double* u = new double[imgRows*imgCols];
    double* v = new double[imgRows*imgCols];
    double* w = new double[imgRows*imgCols];
    double ZMin = *min_element(zIn,zIn+imgRows*imgCols); double ZMax = *max_element(zIn,zIn+imgRows*imgCols);
    double ZDelta = (ZMax - ZMin);
    if(ZDelta == 0.0)
        ZDelta = 1.0;
    double uMin = *min_element(uIn,uIn+imgRows*imgCols); double uMax = *max_element(uIn,uIn+imgRows*imgCols);
    double uDelta = (uMax - uMin);
    if(uDelta == 0.0)
        uDelta = 1.0;
    double vMin = *min_element(vIn,vIn+imgRows*imgCols); double vMax = *max_element(vIn,vIn+imgRows*imgCols);
    double vDelta = (vMax - vMin);
    if(vDelta == 0.0)
        vDelta = 1.0;
    double wMin = *min_element(wIn,wIn+imgRows*imgCols); double wMax = *max_element(wIn,wIn+imgRows*imgCols);
    double wDelta = (wMax - wMin);
    if(wDelta == 0.0)
        wDelta = 1.0;
    for(int idx = 0; idx < imgRows*imgCols; ++idx)
    {
        Z[idx] = (zIn[idx] - ZMin)/ZDelta;
        u[idx] = (uIn[idx] - uMin)/uDelta;
        v[idx] = (vIn[idx] - vMin)/vDelta;
        w[idx] = (wIn[idx] - wMin)/wDelta;
    }

    agglomerativeClustering(imgRows, imgCols, imgChannels, imgLabels, segmentation, img, Z, u, v, w, 0, true, mergeTh);

    delete [] img; delete [] Z; delete [] u; delete [] v; delete [] w;
}



