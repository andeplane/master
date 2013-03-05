#pragma once
#include <armadillo>
using namespace arma;

class Image
{
public:
    Image();
    mat img;
    mat loadBMP(char* filename);
    mat readBMP(const char *filename);
};
