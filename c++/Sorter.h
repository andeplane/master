#ifndef SORTER_H
#define SORTER_H

class System;

#include "System.h"
using namespace arma;

class Sorter {

public:
    // Class data (sorting lists)
    int *cell_n;
    int *index;
    int *Xref;
    
    int ncell;

    System *system;

    // Regular Constructor. 
    Sorter(int ncell, System *system);
    void sort();
};

#endif