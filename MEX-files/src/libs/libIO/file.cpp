#include <iostream>
#include <fstream>

#include "file.hpp"

using namespace std;

bool fexists(const char *filename)
{
    ifstream ifile(filename);
    return ifile;
}
