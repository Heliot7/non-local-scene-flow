#ifndef _SHOWIMAGEMATLAB_
#define _SHOWIMAGEMATLAB_

enum typeMat {CHAR = 1, DOUBLE = 8};

void showImage(unsigned char* img, int imgRows, int imgCols, int nChannels);
void showMatrix(double* mat, int matRows, int matCols);

#endif

