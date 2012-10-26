#ifndef SORTER_H
#define SORTER_H

class System;

#include "system.h"
using namespace arma;

class Sorter {

public:
    System *system;

    // Regular Constructor. 
    Sorter(System *system);
    int sort_system();
};

#endif
