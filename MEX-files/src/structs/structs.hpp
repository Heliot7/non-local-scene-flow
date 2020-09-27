/* Variable wrappers */
struct P2D
{
    double x, y;
    double row, col;
    void initZero() { x = 0.0; y = 0.0; row = 0.0; col = 0.0;}
    void initMinus() { x = -1.0; y = -1.0; row = -1.0; col = -1.0; }
};

struct P3D
{
    double X, Y, Z;
    void initZero() { X = 0.0; Y = 0.0; Z = 0.0; }
};

struct Increments
{
  double dZ, du, dv, dw;
  void initZero() { dZ = 0.0; du = 0.0; dv = 0.0; dw = 0.0;}
};

struct Unknowns
{
  double Z, u, v, w;
  // For clearness when dividing data and smoothness terms.
  double uvw;
  // For clearness when using the konstant part in an equation.
  double konst;
  void initZero() { Z = 0.0; u = 0.0; v = 0.0; w = 0.0; uvw = 0.0; konst = 0.0;}
};

struct ListUnknowns
{
  double* Z;
  double* u;
  double* v;
  double* w;
};

struct BrightnessConstancy
{
    double t0, t1, temp, cross_i1, cross_1i;
    void initZero() { t0 = 0.0; t1 = 0.0; temp = 0.0; cross_i1 = 0.0; cross_1i = 0.0;}
    void initOne() { t0 = 1.0; t1 = 1.0; temp = 1.0; cross_i1 = 1.0; cross_1i = 1.0;}
};

struct MatricesAb
{
    // matrix "a"
    double a[16]; // = a[4][4], but representing in 1D array for MATLAB correctness
    // matrix "b"
    double b[4];
    void initZero()
    {
        for(int i = 0; i < 4; ++i)
        {
            a[i*4] = 0.0;
            a[i*4+1] = 0.0;
            a[i*4+2] = 0.0;
            a[i*4+3] = 0.0;
            b[i] = 0.0;
        }
    }
};

struct DiffNeighbours
{
    double up, down, left, right;
    void initZero() { up = 0.0; down = 0.0; left = 0.0; right = 0.0;}
};

struct Weights
{
    double alpha, alphaZ, alphaUV, alphaW, alphaRig;
    double muZ, muUVW, muS;
    void initZero()
    {
        alpha = 0.0; alphaZ = 0.0; alphaUV = 0.0; alphaW = 0.0; alphaRig = 0.0;
        muZ = 0.0; muUVW = 0.0; muS = 0.0;
    }
};

struct Pixel
{
    int ROW, COL, idx, X;
    void initZero()
    {
        ROW = 0; COL = 0; idx = 0; X = 0;
    }
};

struct pCoef
{
    int ROW, COL, idx, X;
};
