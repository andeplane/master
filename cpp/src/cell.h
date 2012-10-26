#pragma once

class System;
#include <random.h>
#include "molecule.h"
#include "system.h"
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
    double f_sum;
    double density;
    double T;
    int momentum_time_steps;
    vec momentum;
    vec momentum_change;

    int collision_pairs;
    int i,j;
    int particles;
    int particle_capacity;
    unsigned int *particle_indices;


    Cell(System *system);
	void reset();
    void prepare();
    void resize(int n);
    int collide(Random *rnd);
	void resetPressureCalculation();
	void sampleStatistics();
};
