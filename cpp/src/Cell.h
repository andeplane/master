#pragma once

class System;
#include <Random.h>
#include "Molecule.h"
#include "System.h"
#include <armadillo>
#include <CVector.h>

using namespace arma;
using namespace std;

class Cell {
public:
	double volume;

	double vr_max;
    System *system;

	double energy;
    double pressure;
	double density;
    vec momentum;
    vec momentum_change;

    int i,j;

	int particles;
    int particle_capacity;
    unsigned int *particle_indices;


    Cell(System *system);
	void reset();
    void resize(int n);
    int collide(Random *rnd);
	void resetPressureCalculation();
	void sampleStatistics();
};
