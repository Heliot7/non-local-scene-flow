#ifndef _SEGMENTATION_
#define _SEGMENTATION_

struct InputKMeans
{
    int numK, numIterKMeans;
    bool onConnectedComponents;
};

struct Neighbours
{
    int n1, n2;
    double dist;
};

class distComparison
{
public:
    distComparison() {}
    bool operator() (const Neighbours& lhs, const Neighbours&rhs) const
    {
        return (lhs.dist > rhs.dist);
    }
};


void KMeans9D(int* segmentation, int imgRows, int imgCols, unsigned char* image, int imgChannels, double* zIn, double* uIn, double* vIn, double* wIn, InputKMeans input);

void superPixelSegmentation(int imgRows, int imgCols, int imgChannels, int* segmentation, unsigned char* img, double* Z, double* u, double* v, double* w);
int twoPassConnectedComponentAlgorithm(int imgRows, int imgCols, int* segmentation);
void agglomerativeClustering(int imgRows, int imgCols, int imgChannels, int imgLabels, int* segmentation,
                             double* img, double* Z, double* u, double* v, double* w, int thSize, bool similarityAwareness, double mergeTh = 999);

void mergeSegments(int imgRows, int imgCols, int imgChannels, int imgLabels, int* segmentation, unsigned char* img, double* Z, double* u, double* v, double* w, double mergeTh);

#endif

