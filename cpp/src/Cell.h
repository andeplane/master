#pragma once

class System;
#include <Random.h>
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
    int collide(Random *rnd);
	void resetPressureCalculation();
	void sampleStatistics();
};
