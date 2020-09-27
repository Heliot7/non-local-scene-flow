#include <iostream>
#include <algorithm>
#include <math.h>

#include "filters.hpp"

using namespace std;

// Applies median filter to a given image.
void medianFilter(double* medI, int imgRows, int imgCols, double* I, int wSize)
{
    int halfW = (int)floor(wSize/2);
    double w[wSize*wSize];
    for(int COL = halfW; COL < imgCols-halfW; ++COL)
        for(int ROW = halfW; ROW < imgRows-halfW; ++ROW)
        {
            int idx = 0;
            for(int r = ROW-halfW; r <= ROW+halfW; ++r)
                for(int c = COL-halfW; c <= COL+halfW; ++c)
                {
                    w[idx] = I[imgRows*c + r];
                    idx++;
                }
            // Sort them.
            sort(w,w+wSize*wSize);
            // Pick up middle on.
            medI[imgRows*COL + ROW] = w[(int)floor((wSize*wSize)/2)+1];
        }

    // Treat Borders by keeping same value
    for(int COL = halfW; COL < imgCols-halfW; ++COL)
    {
        for(int ROW = 0; ROW < halfW; ++ROW)
            medI[imgRows*COL + ROW] = medI[imgRows*COL + halfW];
        for(int ROW = imgRows-halfW; ROW < imgRows; ++ROW)
            medI[imgRows*COL + ROW] = medI[imgRows*COL + imgRows-halfW-1];
    }

    for(int ROW = 0; ROW < imgRows; ++ROW)
    {
        for(int COL = 0; COL < halfW; ++COL)
            medI[imgRows*COL + ROW] = medI[imgRows*halfW + ROW];
        for(int COL = imgCols-halfW; COL < imgCols; ++COL)
            medI[imgRows*COL + ROW] = medI[imgRows*(imgCols-halfW-1) + ROW];
    }
}
