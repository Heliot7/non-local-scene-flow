#include <iostream>
#include <math.h>

#include "metric.hpp"

double norm(double d1, double d2)
{
    int hola = 1;
    hola++;
    std::cout << "hola" << std::endl;
    return sqrt(d1*d1 + d2*d2) + hola;
}
