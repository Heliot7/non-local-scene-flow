#include <math.h>
#include "matrix.hpp"

// Private function: Computes the inverse of a 3x3 matrix
//  0   3   6
//  1   4   7
//  2   5   8
void inverse3x3(double* M)
{
    // Aux matrix M;
    double a = M[0]; double b = M[3]; double c = M[6];
    double d = M[1]; double e = M[4]; double f = M[7];
    double g = M[2]; double h = M[5]; double i = M[8];

    // Denominator
    double det = det3x3(M);

    // Numerator + Den
    M[0] = (e*i - f*h) / det;
    M[1] = (f*g - d*i) / det;
    M[2] = (d*h - e*g) / det;
    M[3] = (c*h - b*i) / det;
    M[4] = (a*i - c*g) / det;
    M[5] = (b*g - a*h) / det;
    M[6] = (b*f - c*e) / det;
    M[7] = (c*d - a*f) / det;
    M[8] = (a*e - b*d) / det;
}

// Private function: Computes the determinant of a 3x3 matrix
double det3x3(double* M)
{
    return M[0]*M[4]*M[8] + M[3]*M[7]*M[2] + M[1]*M[5]*M[6] - M[6]*M[4]*M[2] - M[3]*M[1]*M[8] - M[7]*M[5]*M[0];
}

double* mulMat(int row1, int col1, double* m1, int row2, int col2, double* m2)
{
    double* resM = new double[row1*col2];

    for(int i = 0; i < row1; ++i)
        for(int j = 0; j < col2; ++j)
        {
            resM[j*row1 + i] = 0.0;
            for(int k = 0; k < col1; ++k)
            {
                resM[j*row1 + i] += m1[k*row1 + i]*m2[j*row2 + k];
            }
        }

    return resM;
}

double* transMat(int row, int col, double* mat)
{
    double* transM = new double[row*col];

    for(int i = 0; i < row; ++i)
        for(int j = 0; j < col; ++j)
        {
            transM[i*col + j] = mat[j*row + i];
        }

    return transM;
}

double* invMat(int n, double* M)
{
    double* invM = new double[n*n];
    double* auxM = new double[n*n];

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            invM[j*n + i] = 0.0;
            auxM[j*n + i] = M[j*n + i];
        }
        invM[i*n + i] = 1.0;
    }

    for(int i = 0; i < n; ++i)
    {
        double temp = auxM[i*n + i];
        if(temp < 0.0)
            temp *= -1;
        int p = i;
        for(int j = i+1; j < n; ++j)
        {
            double tem = 0.0;
            if(auxM[i*n + j] < 0.0)
                tem = -auxM[i*n + j];
            else
                tem = auxM[i*n + j];
            if(temp < 0.0)
                temp *= -1;
            if(tem > temp)
            {
                p = j;
                temp = auxM[i*n + j];
            }

        }
        // Row exchange in both the matrix.
        for(int j = 0; j < n; ++j)
        {
            // Swap 1
            double temp1 = auxM[j*n + i];
            auxM[j*n + i] = auxM[j*n + p];
            auxM[j*n + p] = temp1;
            // Swap 2
            double temp2 = invM[j*n + i];
            invM[j*n + i] = invM[j*n + p];
            invM[j*n + p] = temp2;
        }
        // Dividing the row by a[i*n + i].
        double temp4 = auxM[i*n + i];
        for(int j = 0; j < n; ++j)
        {
            auxM[j*n + i] /= temp4;
            invM[j*n + i] /= temp4;
        }
        // Making other elemes 0 to make 'auxM' the Id and 'invM' the INV.
        for(int q = 0; q < n; ++q)
        {
            if(q == i)
                continue;
            double temp5 = auxM[i*n + q];
            for(int j = 0; j < n; ++j)
            {
                auxM[j*n + q] = auxM[j*n + q] - (temp5*auxM[j*n + i]);
                invM[j*n + q] = invM[j*n + q] - (temp5*invM[j*n + i]);
            }
        }

    }

    delete[] auxM;
    return invM;
}


