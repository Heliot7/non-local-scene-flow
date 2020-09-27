#include <math.h>
#include "differentiation.hpp"

using namespace std;

// Private function: Gradient in X and Y dimensions from given images (finite 1st order central differences).
void gradient(int imgRows, int imgCols, unsigned char* image, double* img_xy)
{
    // Inner pixels
    for(int COL = 1; COL < imgCols-1; ++COL)
        for(int ROW = 1; ROW < imgRows-1; ++ROW)
        {
            // Gradient_x
            img_xy[imgRows*COL + ROW] = (image[imgRows*(COL+1) + ROW] - image[imgRows*(COL-1) + ROW])/2.0;
            // Gradient_y
            img_xy[imgRows*COL + ROW + imgRows*imgCols] = (image[imgRows*COL + ROW+1] - image[imgRows*COL + ROW-1])/2.0;
        }

    // Treat borders by copying adjacent gradients.
    for(int COL = 1; COL < imgCols-1; ++COL)
    {
        // Top
        img_xy[imgRows*COL] = img_xy[imgRows*COL + 1];
        img_xy[imgRows*COL + imgRows*imgCols] = img_xy[imgRows*COL + 1 + imgRows*imgCols];
        // Bottom
        img_xy[imgRows*COL + imgRows-1] = img_xy[imgRows*COL + imgRows-2];
        img_xy[imgRows*COL + imgRows-1 + imgRows*imgCols] = img_xy[imgRows*COL + imgRows-2 + imgRows*imgCols];
    }
    for(int ROW = 0; ROW < imgRows; ++ROW)
    {
        // Left
        img_xy[ROW] = img_xy[imgRows + ROW];
        img_xy[ROW + imgRows*imgCols] = img_xy[imgRows + ROW + imgRows*imgCols];
        // Right
        img_xy[imgRows*(imgCols-1) + ROW] = img_xy[imgRows*(imgCols-2) + ROW];
        img_xy[imgRows*(imgCols-1) + ROW + imgRows*imgCols] = img_xy[imgRows*(imgCols-2) + ROW + imgRows*imgCols];
    }
}

void gradientDouble(int imgRows, int imgCols, double* matrix, double* dX, double* dY)
{
    // Inner pixels
    for(int COL = 1; COL < imgCols-1; ++COL)
        for(int ROW = 1; ROW < imgRows-1; ++ROW)
        {
            // Gradient_x
            dX[imgRows*COL + ROW] = (matrix[imgRows*(COL+1) + ROW] - matrix[imgRows*(COL-1) + ROW])/2.0;
            // Gradient_y
            dY[imgRows*COL + ROW] = (matrix[imgRows*COL + ROW+1] - matrix[imgRows*COL + ROW-1])/2.0;
        }

    // Treat borders by copying adjacent gradients.
    for(int COL = 1; COL < imgCols-1; ++COL)
    {
        // Top
        dX[imgRows*COL] = dX[imgRows*COL + 1];
        dY[imgRows*COL] = dY[imgRows*COL + 1];
        // Bottom
        dX[imgRows*COL + imgRows-1] = dX[imgRows*COL + imgRows-2];
        dY[imgRows*COL + imgRows-1] = dY[imgRows*COL + imgRows-2];
    }
    for(int ROW = 0; ROW < imgRows; ++ROW)
    {
        // Left
        dX[ROW] = dX[imgRows + ROW];
        dY[ROW] = dY[imgRows + ROW];
        // Right
        dX[imgRows*(imgCols-1) + ROW] = dX[imgRows*(imgCols-2) + ROW];
        dY[imgRows*(imgCols-1) + ROW] = dY[imgRows*(imgCols-2) + ROW];
    }
}
