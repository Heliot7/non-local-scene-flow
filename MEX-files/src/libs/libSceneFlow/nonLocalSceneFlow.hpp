#ifndef _NONLOCALSCENEFLOW_
#define _NONLOCALSCENEFLOW_

#include "structs.hpp"
#include <iostream>
#include <string>

struct Input
{
    // Global Paremeters:
    bool onSequence;
    std::string dataset;
    std::string imgExt;
    int numCameras;
    bool onStore;
    bool onDisplay;
    bool onResidual;

    // 3D Scene Flow Parameters:
    // - Multiresolution pyramid.
    int numLevels;              // Number of levesl
    int initLevel;              // Coarsest level
    int endLevel;               // Finest level
    double factor;              // Downsampling factor (0.93)
    // - 2 Fixed point iterations.
    int maxOut;                 // Maximum #loops in outer iteration (2)
    int maxIn;                  // Maximum #loops in inner iteration (15)
    int epsOut;                 // epsilon (cutoff value) to exit outer loops (0.01)
    int epsIn;                  // epsilon (cutoff value) to exit outer loops (0.01)
    // - Initialisation of unknowns.
    int Zinit;                  // Z initialisation.
    bool onZinit;               // Disparity algorithim inits Z or not.
    // - Occlusions parameters.
    bool onOcclusions;          // Occlusion treatment is activated or not
    int border;                 // Number of pixels from the border occluded
    bool onSmoothing;           // Apply gaussian smooth to warped images
    // - Weights for divergence coefficients.
    double alphaZ;
    double alphaUV;
    double alphaW;
    double mu;
    double muZ;
    double muUVW;
    // - iterative solver parameters.
    std::string itSolver;       // Type of iterative solver
    std::string preconditioner; // Type of preconditioner (SSOR, ILU, ILUT, ILUM)
    int ILUp;                   // Level of fill-in in ILU preconditioners
    int maxIter;                // Stablish max. #iterations
    double weightSOR;           // Weight factor
    double errorSolver;         // Maximum relative residual
    // - Epsilon Psi robust functions
    double epsSmooth;           // Epsilon for diffusivity
    double epsData;             // Epsilon for data robust function

    // Segmentation Based:
    // - Segmentation:
    bool onSegmentationBased;
    double proportionKmax;
    double distanceKmax;
    int numBins;
    int numIterKMeans;
    int maxTry;
    // - Disparity:
    double maxDisparity;     // 20-25% of pixel width
    // - Plane fitting:
    int maxIterPlane;
    // - Penalties among segments:
    double pImpact;
    // - Window size: must be half size -> (w*2+1) to ensure odd number.
    double aggregatedWindowSize;
    double dispLambda;
};

struct CameraParameters
{
    // Intrinsic parameters
    double C[9];
    double focalLength[2];
    double principalPoint[2];
    // Extrinsic parameters
    double R[9];
    double T[3];
    double COP[3];
};

struct Result
{

};

class NonLocalSceneFlow
{
    public:

        // Constructor
        NonLocalSceneFlow(Input pInput): input(pInput) {}
        // Triggers the non-local scene flow method
        void run();

    private:

        // Input parameters
        Input input;
        // Unknowns
        double* Z, u, v, w;
};

// Wraps the run method within a MexFile
extern Result nonLocalSceneFlowMatlab(Input input);

// getXYfromZ
extern void getXYfromZ(int imgRows, int imgCols, double* X, double* Y, double* Z, double* focals, double* principals, double* COP);

// projection3Dto2D
extern P2D projection3Dto2D(double* projMatrix, double* P);

// buildProjectionMaps
extern void buildProjectionMaps(double* projMap, int numCameras, int imgRows, int imgCols, double** projMatrices, double* X, double* Y, double* Z, double* u = NULL, double* v = NULL, double* w = NULL);
P2D projection3Dto2D(double* projMatrix, double X, double Y, double Z);

// buildOcclusionMaps
void buildOcclusionMaps(bool* occlusionMap, int imgRows, int imgCols, int numCameras, double* projImg, double* X, double* Y, double* Z, double* u, double* v, double* w, double** COP);
bool isPxlOutOfBounds(int imgRows, int imgCols, P2D pxl);

// buildOcclusionBinary
extern void buildOcclusionBinary(bool* occBin, int numCameras, int imgRows, int imgCols, double* t0ProjImg, double* t1ProjImg, bool* t0occlusionMaps, bool* t1occlusionMaps);

// buildWarpedImage
extern void buildWarpedImage(double* warpedImgs, int numCameras, int imgRows, int imgCols, int imgChannels, unsigned char** images, double* projImgs, bool* occlusionMaps);
double bilinearInterpolation(double ROW, double COL, P2D p11, double i11, P2D p12, double i12, P2D p21, double i21, P2D p22, double i22);

// computePartialUnknowns
extern void computePartialUnknowns(double* partialUnknowns, int numCameras, int imgRows, int imgCols, int imgChannels, double** projMatrices, double* focals, double* principals, double* camPos,
                                           double* Z, double* Z_x, double* Z_y, double* u, double* u_x, double* u_y, double* v, double* v_x, double* v_y, double* w, double* w_x, double* w_y,
                                           double* warpedImgs_x, double* warpedImgs_y, int numFrame);
void buildJacobianT0(double* J, double* M, P2D f, P2D cam, double X, double Y, double Z, double Z_x, double Z_y);
void buildJacobianT1(double* J, double* M, P2D f, P2D cam, double X, double Y, double Z, double u, double v, double w, double Z_x, double Z_y, double u_x, double u_y, double v_x, double v_y, double w_x, double w_y);
P2D p_iWarpedGradient(P2D gradPxlRefImg, double* J);
P2D t0PxlPartial_Z(double* M, P2D cam, double X, double Y, double Z);
P2D t1PxlPartial_Z(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w);
P2D t1PxlPartial_u(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w);
P2D t1PxlPartial_v(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w);
P2D t1PxlPartial_w(double* M, P2D cam, double X, double Y, double Z, double u, double v, double w);
P2D getXYnoZ(int pox, int posy, double focalX, double focalY, double principalX, double principalY);
void inverse2x2(double* invM, double* M);

// buildDiffusivitySmooth
extern void buildDiffusivitySmooth(double* diffSmooth, int imgRows, int imgCols, double* Z_x, double* Z_y, double* u_x, double* u_y, double* v_x, double* v_y, double* w_x, double* w_y,
                                        double* dZ_x, double* dZ_y, double* du_x, double* du_y, double* dv_x, double* dv_y, double* dw_x, double* dw_y,
                                        double epsilon, double muZ, double muUVW);

// buildRobustnessData
extern void buildRobustnessData(double* psiData, int numCameras, int imgRows, int imgCols, int imgChannels, double* dZ, double* du, double* dv, double* dw,
                                                 double* t0WarpedImgs, double* t1WarpedImgs, double* t0I_unknown, double* t1I_unknown, double epsilon);
double psiDeriv(double s, double epsilon);

// buildMatricesAb
extern void buildMatricesAb(double* mA, double* mB, int numCameras, int imgRows, int imgCols, int imgChannels, double* Z, double* u, double* v, double* w,
                            double* t0WarpedImgs, double* t1WarpedImgs, double* t0I_unknown, double* t1I_unknown, double* psiData, double* psiSmooth, bool* occlusionBinary,
                            double LEVEL, int numLevels, double factor, Weights weights, double* segmentPenalties);
DiffNeighbours computeDiffusivityZ(int imgRows, int imgCols, double* psiSmooth, int ROW, int COL);
DiffNeighbours computeDiffusivityUVW(int imgRows, int imgCols, double* psiSmooth, int ROW, int COL);

/* Iterative Solvers */
extern void iterativeSolver(double* dZ, double* du, double* dv, double* dw, double* evolResidual, char* itSolver, char* itPrecon, int imgRows, int imgCols, double* a, double* b, double* diffSmooth,
                                   double wSOR, int maxIter, double errorSolver, int LEVEL, int numLevels, double factor, Weights weights, double* segmentPenalties);
Increments treatBorders(int imgRows, int imgCols, int row, int col, double* dZ, double* du, double* dv, double* dw, int limit);

// SOR
extern Increments SOR(int imgRows, int imgCols, double* a, double* b, int ROW, int COL, double* incr, double* div_d, double w);
void residualSystem(int imgRows, int imgCols, double* residual, double* a, int ROW, int COL, double* b, double* x);
// CG & CR
void CG_incr_res(double* incr, double* r, double* p, double* Ap, double alpha);
void CG_dir(double* r, double* d, double beta);
void mulAp(double* Ap, double* A, double* p);

// Preconditioning
void computePreconditioner(double** precA, int imgRows, int imgCols, int limit, char* preconditioner, double** A);
void computeILU(int imgRows, int imgCols, int limit, double** M);
void applyPreconditioner(int imgRows, int imgCols, int limit, double** z, double** r, double** M);

#endif

