#include <math.h>
#include "colour.hpp"

// Converts an image from RGB to greyscale (storage of 1D intensties in the first RGB channel
void rgbToGreyscale(int imgRows, int imgCols, unsigned char* imgGrey, unsigned char* img)
{
    int size = imgRows*imgCols;
    for(int idx = 0; idx < size; ++idx)
        imgGrey[idx] = (unsigned char)(0.21 * img[idx + size*0] + 0.71 * img[idx + size*1] + 0.07 * img[idx + size*2]);
}
