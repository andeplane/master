#ifndef CELL_H
#define CELL_H

class System;

#include "Molecule.h"
#include "System.h"

using namespace std;

class Cell {
private:
	
public:
	double vr_max;
	double selxtra;
	System *system;
	Cell(System *system);

};

#endif