#ifndef CELL_H
#define CELL_H

class System;

#include "Molecule.h"
#include "System.h"
#include <armadillo>

using namespace arma;
using namespace std;

class Cell {
private:
	
public:
	int index;
	int firstParticleIndex;
	double volume;

	double vr_max;
	System *system;

	double energy;
	double density;
	vec momentum;

	int particles;

	Cell(System *system);
	void reset();
	int collide();
	void resetPressureCalculation();
	void sampleStatistics();
};

#endif