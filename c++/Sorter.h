#ifndef SORTER_H
#define SORTER_H

class System;

#include "System.h"
using namespace arma;

class Sorter {

public:
    int *Xref;
    System *system;

    // Regular Constructor. 
    Sorter(System *system);
    void sort();
};

#endif