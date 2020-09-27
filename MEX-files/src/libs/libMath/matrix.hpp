#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>
#include <string>
#include <assert.h>

template <class T>
class Matrix
{
    public:

        // Constructors
        Matrix(): imgRows(0), imgCols(0), imgChannels(0), matrix(NULL) {}
        Matrix(const Matrix<T>& copyMatrix)
        {
            matrix = NULL;
            (*this) = copyMatrix;
        }
        Matrix(int rows, int cols, int channels = 1, bool bInit = false, int iInit = 0): imgRows(rows), imgCols(cols), imgChannels(channels)
        {
            matrix = NULL;
            createMatrix(rows, cols, channels);
            // Init to zero if not iInit specified
            if(bInit)
                init(iInit);
        }

        // Destructor
        ~Matrix() { clean(); }

        // Init
        void init(T initValue)
        {
            unsigned int numElems = imgRows*imgCols*imgChannels;
            for(unsigned int idx = 0; idx < numElems; ++idx)
                matrix[idx] = initValue;
        }

        // Free allocated memory
        void clean()
        {
            if(matrix != NULL)
            {
                delete [] matrix;
                matrix = NULL;
            }
        }

        // Print
        void print(std::string name)
        {
            // Prints only first channels
            std::cout << "Matrix " << name << ":" << std::endl;
            for(int row = 0; row < getRows(); ++row)
            {
                for(int col = 0; col < getCols(); ++col)
                    std::cout << getElem(row,col) << " ";
                std::cout << std::endl;
            }
        }

        // Getters and Setters

        // - Matrix
        void createMatrix(int rows, int cols, int channels = 1) { imgRows = rows; imgCols = cols; imgChannels = channels; clean(); matrix = new T[imgRows*imgCols*imgChannels]; }
        void setMatrix(T* pMatrix) { matrix = pMatrix; }
        T* getMatrix() const { return matrix; }
        bool isNull() const { return (matrix == NULL) ? true : false; }

        // - Elements of the matrix
        T getElem(int row, int col, int channel = 0) const { return matrix[imgRows*col + row + channel*imgRows*imgCols]; }
        void setElem(T elem, int row, int col, int channel = 0) { matrix[imgRows*col + row + channel*imgRows*imgCols] = elem; }

        // - Other class attributes
        template <class U>
        bool isSameSize(Matrix<U>& m) { return (getRows() == m.getRows() && getCols() == m.getCols() && getChannels() == m.getChannels()); }
        int getNumElems() const { return getRows() * getCols() * getChannels(); }
        int getRows() const { return imgRows; }
        void setRows(int pRows) { imgRows = pRows; }
        int getCols() const { return imgCols; }
        void setCols(int pCols) { imgCols = pCols; }
        int getChannels() const { return imgChannels; }
        void setChannels(int pChannels) { imgChannels = pChannels; }

        // Overloading operators

        // - Assignemnt whole matrix (by copy element by element)
        Matrix<T>& operator=(const Matrix<T>& rhs)
        {
            // Same object ?
            if (this == &rhs)
                return *this;

            if(getMatrix() != NULL)
                delete [] getMatrix();
            setRows(rhs.getRows());
            setCols(rhs.getCols());
            setChannels(rhs.getChannels());
            setMatrix(new T[getRows()*getCols()*getChannels()]);

            T* lhsMatrix = this->getMatrix();
            T* rhsMatrix = rhs.getMatrix();

            // Copy elements of the matrix (rhs)
            int numElems = getRows() * getCols() * getChannels();
            for(int idx = 0; idx < numElems; ++idx)
                lhsMatrix[idx] = rhsMatrix[idx];

            return *this;
        }

        template <class U>
        Matrix<T>& operator=(const Matrix<U>& rhs)
        {
            if(getMatrix() != NULL)
                delete [] getMatrix();
            setRows(rhs.getRows());
            setCols(rhs.getCols());
            setChannels(rhs.getChannels());
            setMatrix(new T[getRows()*getCols()*getChannels()]);

            T* lhsMatrix = this->getMatrix();
            U* rhsMatrix = rhs.getMatrix();

            // Copy elements of the matrix (rhs)
            int numElems = getRows() * getCols() * getChannels();
            for(int idx = 0; idx < numElems; ++idx)
                lhsMatrix[idx] = static_cast<T>(rhsMatrix[idx]);

            return *this;
        }

        // - Parenthesis for better element processing
        const T operator()(const int row, const int col, const int channel = 0) const
        {
            return matrix[imgRows*col + row + channel*imgRows*imgCols];
        }
        T& operator()(const int row, const int col, const int channel = 0)
        {
            return matrix[imgRows*col + row + channel*imgRows*imgCols];
        }

        // Arithmetic operations (static)

        // - Addition (unary)
        template <class U>
        static inline void add(Matrix<T>& Mout, Matrix<U>& Min)
        {
            assert(Mout.isSameSize(Min));

            T* matrixOut = Mout.getMatrix();
            U* matrixIn = Min.getMatrix();

            int numElems = Mout.getNumElems();
            for(int idx = 0; idx < numElems; ++idx)
                matrixOut[idx] += static_cast<T>(matrixIn[idx]);
        }

        // - Addition (binary)
        template <class U, class V>
        static inline void add(Matrix<T>& Mout, Matrix<U>& M1, Matrix<V>& M2)
        {
            assert(M1.isSameSize(M2));

            if(!Mout.isNull() && !Mout.isSameSize(M1))
                Mout.createMatrix(M1.getRows(), M1.getCols(), M1.getChannels());

            T* matrixOut = Mout.getMatrix();
            U* matrixIn1 = M1.getMatrix();
            V* matrixIn2 = M2.getMatrix();

            int numElems = Mout.getNumElems();
            for(int idx = 0; idx < numElems; ++idx)
                matrixOut[idx] = static_cast<T>(matrixIn1[idx] + matrixIn2[idx]);
        }

        // - Substraction (unary)
        template <class U>
        static inline void sub(Matrix<T>& Mout, Matrix<U>& Min)
        {
            assert(Mout.isSameSize(Min));

            T* matrixOut = Mout.getMatrix();
            U* matrixIn = Min.getMatrix();

            int numElems = Mout.getNumElems();
            for(int idx = 0; idx < numElems; ++idx)
                matrixOut[idx] -= static_cast<T>(matrixIn[idx]);
        }

        // - Substraction (binary)
        template <class U, class V>
        static inline void sub(Matrix<T>& Mout, Matrix<U>& M1, Matrix<V>& M2)
        {
            assert(M1.isSameSize(M2));

            if(!Mout.isNull() && !Mout.isSameSize(M1))
                Mout.createMatrix(M1.getRows(), M1.getCols(), M1.getChannels());

            T* matrixOut = Mout.getMatrix();
            U* matrixIn1 = M1.getMatrix();
            V* matrixIn2 = M2.getMatrix();

            int numElems = Mout.getNumElems();
            for(int idx = 0; idx < numElems; ++idx)
                matrixOut[idx] = static_cast<T>(matrixIn1[idx] - matrixIn2[idx]);
        }

        // - Multiplication (binary)
        template <class U, class V>
        static inline void mul(Matrix<T>& Mout, Matrix<U>& M1, Matrix<V>& M2)
        {
            int numRows1 = M1.getRows();
            int numRows2 = M2.getRows();
            int numCols1 = M1.getCols();
            int numCols2 = M2.getCols();
            int numChannels = M1.getChannels();

            assert(numCols1 == numRows2 && M1.getChannels() == M2.getChannels());

            Mout.createMatrix(numRows1, numCols2, numChannels);
            for(int chn = 0; chn < numChannels; ++chn)
            {
                for(int row = 0; row < numRows1; ++row)
                    for(int col = 0; col < numCols2; ++col)
                    {
                        Mout(row,col,chn) = 0;
                        for(int colLine = 0; colLine < numCols1; ++colLine)
                            Mout(row,col,chn) += static_cast<T>(M1(row,colLine,chn) * M2(colLine,col,chn));
                    }
            }
        }

        // - Transpose
        template <class U>
        static inline void transpose(Matrix<T>& Mout, Matrix<U>& Min)
        {
            if( !(Mout.getRows() == Min.getCols() && Mout.getCols() == Min.getRows()) )
                Mout.createMatrix(Min.getCols(), Min.getRows(), Min.getChannels());

            int numRows = Mout.getRows();
            int numCols = Mout.getCols();
            int numChannels = Mout.getChannels();
            for(int chn = 0; chn < numChannels; ++chn)
                for(int row = 0; row < numRows; ++row)
                    for(int col = 0; col < numCols; ++col)
                        Mout(row,col,chn) = static_cast<T>( Min(col,row,chn) );
        }

        // - Determinant (2x2)
        T det2()
        {
            assert(getRows() == 2 && getCols() == 2);
            Matrix<T>& mat = (*this);
            return mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0);
        }

        // - Determinant (3x3)
        T det3()
        {
            assert(getRows() == 3 && getCols() == 3);
            Matrix<T>& mat = (*this);
            return mat(0,0) * mat(1,1) * mat(2,2) + mat(0,1) * mat(1,2) * mat(2,0) + mat(0,2) * mat(1,0) * mat(2,1)
                    - mat(0,2) * mat(1,1) * mat(2,0) - mat(0,1) * mat(1,0) * mat(2,2) - mat(0,0) * mat(1,2) * mat(2,1);
        }

        // - Inverse (2x2)
        template <class U>
        static inline void inverse2(Matrix<T>& Mout, Matrix<U>& Min)
        {
            assert(Min.getRows() == Min.getCols());

            Mout.createMatrix(Min.getRows(), Min.getCols());

            T det = Min.det2();

            Mout(0,0) = static_cast<T>( Min(1,1)) / det;
            Mout(0,1) = static_cast<T>(-Min(0,1)) / det;
            Mout(1,0) = static_cast<T>(-Min(1,0)) / det;
            Mout(1,1) = static_cast<T>( Min(0,0)) / det;
        }

        // - Inverse (3x3)
        template <class U>
        static inline void inverse3(Matrix<T>& Mout, Matrix<U>& Min)
        {
            assert(Min.getRows() == Min.getCols());

            Mout.createMatrix(Min.getRows(), Min.getCols());

            T det = Min.det3();

            Mout(0,0) = ( Min(1,1)*Min(2,2) - Min(1,2)*Min(2,1) ) / det;
            Mout(0,1) = ( Min(0,2)*Min(2,1) - Min(0,1)*Min(2,2) ) / det;
            Mout(0,2) = ( Min(0,1)*Min(1,2) - Min(0,2)*Min(1,1) ) / det;
            Mout(1,0) = ( Min(1,2)*Min(2,0) - Min(1,0)*Min(2,2) ) / det;
            Mout(1,1) = ( Min(0,0)*Min(2,2) - Min(0,2)*Min(2,0) ) / det;
            Mout(1,2) = ( Min(0,2)*Min(1,0) - Min(0,0)*Min(1,2) ) / det;
            Mout(2,0) = ( Min(1,0)*Min(2,1) - Min(1,1)*Min(2,0) ) / det;
            Mout(2,1) = ( Min(2,0)*Min(0,1) - Min(0,0)*Min(2,1) ) / det;
            Mout(2,2) = ( Min(0,0)*Min(1,1) - Min(0,1)*Min(1,0) ) / det;
        }

    private:

        T* matrix;
        int imgRows, imgCols, imgChannels;
};

// Procedural functions:
void inverse3x3(double* M);
double det3x3(double* M);

double* mulMat(int row1, int col1, double* m1, int row2, int col2, double* m2);
double* transMat(int row, int col, double* mat);
double* invMat(int n, double* M);

#endif

