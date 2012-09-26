#ifndef CELL_H
#define CELL_H

class System;

#include "Molecule.h"
#include "System.h"

using namespace std;

class Cell {
private:
	
public:
	int index;
	int firstParticleIndex;

	double vr_max;
	double selxtra;
	System *system;

	double timeForPressureReset;
	double delta_v;

	int particlesInCell;

	Cell(System *system);
	void reset();
	int collide();
	void resetPressureCalculation();
};

#endif